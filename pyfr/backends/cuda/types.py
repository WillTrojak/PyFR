# -*- coding: utf-8 -*-

import pycuda.driver as cuda
import numpy as np

import pyfr.backends.base as base

from pyfr.backends.cuda.util import memcpy2d_htod, memcpy2d_dtoh


class CudaMatrixBase(base.MatrixBase):
    def __init__(self, backend, dtype, ioshape, initval, iopacking, tags):
        super(CudaMatrixBase, self).__init__(backend, ioshape, iopacking, tags)

        # Data type info
        self.dtype = dtype
        self.itemsize = np.dtype(dtype).itemsize

        # Compute the size, in bytes, of the minor dimension
        colsz = self.ncol*self.itemsize

        if 'nopad' not in tags:
            # Allocate a 2D array aligned to the major dimension
            self.data, self.pitch = cuda.mem_alloc_pitch(colsz, self.nrow,
                                                         self.itemsize)
            self._nbytes = self.nrow*self.pitch

            # Ensure that the pitch is a multiple of itemsize
            assert (self.pitch % self.itemsize) == 0
        else:
            # Allocate a standard, tighly packed, array
            self._nbytes = colsz*self.nrow
            self.data = cuda.mem_alloc(self._nbytes)
            self.pitch = colsz

        self.leaddim = self.pitch / self.itemsize
        self.traits = (self.nrow, self.leaddim, self.dtype)

        # Process any initial values
        if initval is not None:
            self._set(self._pack(initval))

    def _get(self):
        # Allocate an empty buffer
        buf = np.empty((self.nrow, self.ncol), dtype=self.dtype)

        # Copy
        memcpy2d_dtoh(buf, self.data, self.pitch, self.ncol*self.itemsize,
                      self.ncol*self.itemsize, self.nrow)

        return buf

    def _set(self, ary):
        nary = np.asanyarray(ary, dtype=self.dtype, order='C')

        # Copy
        memcpy2d_htod(self.data, nary, self.ncol*self.itemsize,
                      self.pitch, self.ncol*self.itemsize, self.nrow)

    def offsetof(self, i, j):
        if i < 0 or i >= self.nrow or j < 0 or j >= self.ncol:
            raise ValueError('Index ({},{}) out of bounds ({},{}))'
                             .format(i, j, self.nrow, self.ncol))

        return self.pitch*i + j*self.itemsize

    def addrof(self, i, j):
        return np.intp(long(self.data) + self.offsetof(i, j))

    @property
    def _as_parameter_(self):
        return long(self.data)

    def __long__(self):
        return long(self.data)

    @property
    def nbytes(self):
        return self._nbytes


class CudaMatrix(CudaMatrixBase, base.Matrix):
    def __init__(self, backend, ioshape, initval, iopacking, tags):
        super(CudaMatrix, self).__init__(backend, backend.fpdtype, ioshape,
                                         initval, iopacking, tags)


class CudaMatrixRSlice(base.MatrixRSlice):
    def __init__(self, backend, mat, p, q):
        super(CudaMatrixRSlice, self).__init__(backend, mat, p, q)

        # Copy over common attributes
        self.dtype, self.itemsize = mat.dtype, mat.itemsize
        self.pitch, self.leaddim = mat.pitch, mat.leaddim

        # Traits are those of a thinner matrix
        self.traits = (self.nrow, self.leaddim, self.dtype)

        # Starting offset of our row
        self._soffset = mat.offsetof(p, 0)

    @property
    def _as_parameter_(self):
        return long(self.parent) + self._soffset

    @property
    def __long__(self):
        return long(self.parent) + self._soffset


class CudaMatrixBank(base.MatrixBank):
    def __init__(self, backend, mats, initbank, tags):
        for m in mats[1:]:
            if m.traits != mats[0].traits:
                raise ValueError('Matrices in a bank must be homogeneous')

        super(CudaMatrixBank, self).__init__(mats, initbank, tags)

    @property
    def _as_parameter_(self):
        return long(self._curr_mat)

    def __long__(self):
        return long(self._curr_mat)


class CudaConstMatrix(CudaMatrixBase, base.ConstMatrix):
    def __init__(self, backend, initval, iopacking, tags):
        ioshape = initval.shape
        super(CudaConstMatrix, self).__init__(backend, backend.fpdtype,
                                              ioshape, initval, iopacking,
                                              tags)


class CudaSparseMatrix(object):
    def __init__(self, backend, initval, tags):
        raise NotImplementedError('SparseMatrix todo!')


class CudaView(base.View):
    def __init__(self, backend, matmap, rcmap, stridemap, vlen, tags):
        self.nrow = nrow = matmap.shape[0]
        self.ncol = ncol = matmap.shape[1]
        self.vlen = vlen

        # For vector views a stridemap is required
        if vlen != 1 and np.any(stridemap == 0):
            raise ValueError('Vector views require a non-zero stride map')

        # Check all of the shapes match up
        if matmap.shape != rcmap.shape[:2] or\
           matmap.shape != stridemap.shape:
            raise TypeError('Invalid matrix shapes')

        # Get the different matrices which we map onto
        self._mats = list(np.unique(matmap))

        # Extract the data type and item size from the first matrix
        self.refdtype = self._mats[0].dtype
        self.refitemsize = self._mats[0].itemsize

        # Validate the matrices
        for m in self._mats:
            if not isinstance(m, (CudaMatrix, CudaConstMatrix)):
                raise TypeError('Incompatible matrix type for view')

            if m.dtype != self.refdtype:
                raise TypeError('Mixed view matrix types are not supported')

        # Go from matrices and row/column indices to addresses
        r, c = rcmap[...,0], rcmap[...,1]
        ptrmap = np.vectorize(lambda m, r, c: m.addrof(r, c))(matmap, r, c)

        self.mapping = CudaMatrixBase(backend, np.intp, (nrow, ncol), ptrmap,
                                      'AoS', tags)
        self.strides = CudaMatrixBase(backend, np.int32, (nrow, ncol),
                                      stridemap, 'AoS',tags)

    @property
    def nbytes(self):
        return self.mapping.nbytes + self.strides.nbytes


class CudaMPIMatrix(CudaMatrix, base.MPIMatrix):
    def __init__(self, backend, ioshape, initval, iopacking, tags):
        # Ensure that our CUDA buffer will not be padded
        ntags = tags | {'nopad'}

        # Call the standard matrix constructor
        super(CudaMPIMatrix, self).__init__(backend, ioshape, initval,
                                            iopacking, ntags)

        # Allocate a page-locked buffer on the host for MPI to send/recv from
        self.hdata = cuda.pagelocked_empty((self.nrow, self.ncol),
                                           self.dtype, 'C')


class CudaMPIView(base.MPIView):
    def __init__(self, backend, matmap, rcmap, stridemap, vlen, tags):
        self.nrow = nrow = matmap.shape[0]
        self.ncol = ncol = matmap.shape[1]
        self.vlen = vlen

        # Create a normal CUDA view
        self.view = backend._view(matmap, rcmap, stridemap, vlen, tags)

        # Now create an MPI matrix so that the view contents may be packed
        self.mpimat = backend._mpi_matrix((nrow, ncol, vlen), None, 'AoS',
                                          tags=tags)

    @property
    def nbytes(self):
        return self.view.nbytes + self.mpimat.nbytes
