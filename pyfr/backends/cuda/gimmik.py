# -*- coding: utf-8 -*-

from ctypes import create_string_buffer
from gimmik import generate_mm, generate_tfmm, generate_tfmm_lines, generate_tfmm_managed
import numpy as np
from math import ceil, floor
import subprocess

from pyfr.backends.base import ComputeKernel, NotSuitableError
from pyfr.backends.cuda.provider import (CUDAKernelProvider,
                                         get_grid_for_block)


class CUDAGiMMiKKernels(CUDAKernelProvider):
    def __init__(self, backend):
        super().__init__(backend)

    def mul(self, a, b, out, alpha=1.0, beta=0.0):
        # Ensure the matrices are compatible
        if a.nrow != out.nrow or a.ncol != b.nrow or b.ncol != out.ncol:
            raise ValueError('Incompatible matrices for out = a*b')

        # Check that A is constant
        if 'const' not in a.tags:
            raise NotSuitableError('GiMMiK requires a constant a matrix')

        # Fetch the matrix and tally up the number of non-zeros
        arr = a.get()
        nnz, nuq = np.count_nonzero(arr), len(np.unique(np.abs(arr)))

        # Check that A is suitable
        if nuq > 28 and nnz / arr.size > 0.15:
            raise NotSuitableError('Matrix is inappropriate for GiMMiK')

        # Generate
        src = generate_mm(arr, dtype=a.dtype, platform='cuda',
                          alpha=alpha, beta=beta)

        # Build
        fun = self._build_kernel('gimmik_mm', src,
                                 [np.int32, np.intp]*2 + [np.int32])
        fun.set_cache_pref(prefer_l1=True)


        # Determine the grid/block
        block = (128, 1, 1)
        grid = get_grid_for_block(block, b.ncol)

        class MulKernel(ComputeKernel):
            def run(self, queue):
                fun.exec_async(grid, block, queue.stream_comp, b.ncol, b,
                               b.leaddim, out, out.leaddim)

        return MulKernel()

    def mul_tensor(self, a, b, out):
        # Ensure the matrices are compatible
        print(f'{a.nrow=} {a.ncol=} {b.nrow=} {b.ncol} {out.nrow=} {out.ncol=}')

        # Check that A is constant
        if 'const' not in a.tags:
            raise NotSuitableError('GiMMiK requires a constant a matrix')

        # Fetch the matrix and tally up the number of non-zeros
        arr = a.get()
        nnz, nuq = np.count_nonzero(arr), len(np.unique(np.abs(arr)))

        (p, k) = np.shape(arr)
        nele = int(b.ncol/13)

        shr_size = 0

        # Generate
        if a.dtype == np.float32:
            if p == 5:
                l1_pref  = False
                shr_pref = True
                src, shr_size = generate_tfmm_lines(arr, ndims=3, nvars=13,
                                    dtype=a.dtype, soasz=32, shr_max=96000,
                                    platform='cuda_tensor_lines_24')
                block = (200, 1, 1)
                grid = (ceil(nele*p*p/block[0]), 1, 1)
            elif p == 4:
                l1_pref  = False
                shr_pref = True
                opargs = {'ld_opt': True, 'st_opt': False, 'intl_opt': False,
                          'pipe_opt': False, 'max_coag': 0, 'compute_size': 0,
                          'shr_op_order': 'grs', 'shr_bdc': True}
                src, shr_size = generate_tfmm_managed(arr, ndims=3, nvars=13, 
                                    soasz=32, dtype=a.dtype, opargs=opargs,
                                    platform='cuda_tfmm_managed', flux='hyperns',
                                    shared_max=67584, ufc_size=4*(2**15), 
                                    block_dim=128)
                block = (128, 1, 1)
                warp_size = 32
                elem_warp = int(warp_size/p)
                elem_block = int(block[0]/warp_size)*elem_warp
                grid = (ceil(nele/elem_block), 1, 1)
            elif p <= 3:
                l1_pref  = False
                shr_pref = True

                src, shr_size = generate_tfmm(arr, ndims=3, nvars=13, dtype=a.dtype,
                                    block_dim=128, soasz=32, flux='hyper',
                                    platform='cuda_tensor_flux12')
                block = (128, 1, 1)
                warp_size = 32
                elem_warp = int(warp_size/p)
                elem_block = int(block[0]/warp_size)*elem_warp
                grid = (ceil(nele/elem_block), 1, 1)        
        elif a.dtype == np.float64:
            if p == 5:
                l1_pref  = False
                shr_pref = True
                src, shr_size = generate_tfmm_lines(arr, ndims=3, nvars=13,
                                    dtype=a.dtype, soasz=32, shr_max=96000,
                                    platform='cuda_tensor_lines_12')
                block = (200, 1, 1)
                grid = (ceil(nele*p*p/block[0]), 1, 1)
            if p == 4:
                l1_pref  = False
                shr_pref = True
                src, shr_size = generate_tfmm_lines(arr, ndims=3, nvars=13,
                                    dtype=a.dtype, soasz=32, shr_max=92160,
                                    platform='cuda_tensor_lines_15')
                block = (192, 1, 1)
                grid = (ceil(nele*p*p/block[0]), 1, 1)
            if p == 3:
                l1_pref  = False
                shr_pref = True

                src, shr_size = generate_tfmm(arr, ndims=3, nvars=13, dtype=a.dtype,
                                    block_dim=128, soasz=32, flux='hyper',
                                    platform='cuda_tensor_flux12')
                block = (128, 1, 1)
                warp_size = 32
                elem_warp = int(warp_size/p)
                elem_block = int(block[0]/warp_size)*elem_warp
                grid = (ceil(nele/elem_block), 1, 1)
                # l1_pref  = False
                # shr_pref = True
                # shr_max = 60416
                # opargs = {'ld_opt': True, 'st_opt': False, 'intl_opt': False,
                #         'pipe_opt': False, 'max_coag': 0, 'compute_size': 0,
                #         'shr_op_order': 'grs', 'shr_bdc': True}
                # src, shr_size = generate_tfmm_managed(arr, ndims=3, nvars=13, 
                #                     soasz=32, dtype=a.dtype, opargs=opargs,
                #                     platform='cuda_tfmm_managed', flux='hyperns',
                #                     shared_max=shr_max, ufc_size=4*(2**15), 
                #                     block_dim=128)
                # block = (128, 1, 1)
                # grid = (ceil(nele*p/block[0]), 1, 1)
            if p == 2:
                l1_pref  = False
                shr_pref = True
                shr_max = 43008
                opargs = {'ld_opt': True, 'st_opt': False, 'intl_opt': False,
                        'pipe_opt': False, 'max_coag': 0, 'compute_size': 0,
                        'shr_op_order': 'grs', 'shr_bdc': True}
                src, shr_size = generate_tfmm_managed(arr, ndims=3, nvars=13, 
                                    soasz=32, dtype=a.dtype, opargs=opargs,
                                    platform='cuda_tfmm_managed', flux='hyperns',
                                    shared_max=shr_max, ufc_size=4*(2**15), 
                                    block_dim=128)
                block = (128, 1, 1)
                grid = (ceil(nele*p/block[0]), 1, 1)

        print(f'{nele=}, {shr_size=}, {p=}, {block=}, {grid=}')

        base_name =  f'tmul_p{p-1}'
        src_name = base_name + '.cu'
        ptx_name = base_name + '.ptx'
        f = open(src_name, 'w')
        f.write(src)
        f.close()

        # Build
        fun = self._build_kernel('gimmik_tfmm', src,
                                 [np.int32, np.intp]*2 + [np.int32])#, ptx=ptx)
        fun.set_cache_pref(prefer_l1=l1_pref, prefer_shared=shr_pref)
        fun.set_shared_size(dynm_shared=shr_size)

        class MulKernel(ComputeKernel):
            def run(self, queue):
                fun.exec_async(grid, block, queue.stream_comp, nele, b,
                               b.leaddim, out, out.leaddim,
                               dynm_shared=shr_size)

        return MulKernel()


    def mul_tensor_source(self, a, b, out):
        # Check that A is constant
        if 'const' not in a.tags:
            raise NotSuitableError('GiMMiK requires a constant a matrix')

        # Fetch the matrix and tally up the number of non-zeros
        arr = a.get()
        nnz, nuq = np.count_nonzero(arr), len(np.unique(np.abs(arr)))

        (p, k) = np.shape(arr)
        nele = int(b.ncol/13)

        shr_size = 0

        # Generate
        if a.dtype == np.float32:
            if p == 5:
                l1_pref  = False
                shr_pref = True
                src, shr_size = generate_tfmm_lines(arr, ndims=3, nvars=13,
                                    dtype=a.dtype, soasz=32, shr_max=96000,
                                    platform='cuda_S_tensor_lines_24')
                block = (200, 1, 1)
                grid = (ceil(nele*p*p/block[0]), 1, 1)
            elif p == 4:
                l1_pref  = False
                shr_pref = True
                opargs = {'ld_opt': True, 'st_opt': False, 'intl_opt': False,
                        'pipe_opt': False, 'max_coag': 0, 'compute_size': 0,
                        'shr_op_order': 'grs', 'shr_bdc': True}
                src, shr_size = generate_tfmm_managed(arr, ndims=3, nvars=13, 
                                    soasz=32, dtype=a.dtype, opargs=opargs,
                                    platform='cuda_tfmm_managed', flux='hyperns',
                                    shared_max=67584, ufc_size=4*(2**15), 
                                    block_dim=128, source_term=True)
                block = (128, 1, 1)
                warp_size = 32
                elem_warp = int(warp_size/p)
                elem_block = int(block[0]/warp_size)*elem_warp
                grid = (ceil(nele/elem_block), 1, 1)
            elif p <= 3:
                l1_pref  = False
                shr_pref = True

                src, shr_size = generate_tfmm(arr, ndims=3, nvars=13, dtype=a.dtype,
                                    block_dim=128, soasz=32, flux='hyper',
                                    platform='cuda_S_tensor_flux12')
                block = (128, 1, 1)
                warp_size = 32
                elem_warp = int(warp_size/p)
                elem_block = int(block[0]/warp_size)*elem_warp
                grid = (ceil(nele/elem_block), 1, 1)
        if a.dtype == np.float64:
            if p == 5:
                l1_pref  = False
                shr_pref = True
                src, shr_size = generate_tfmm_lines(arr, ndims=3, nvars=13,
                                    dtype=a.dtype, soasz=32, shr_max=96000,
                                    platform='cuda_S_tensor_lines_12')
                block = (200, 1, 1)
                grid = (ceil(nele*p*p/block[0]), 1, 1)
            if p == 4:
                l1_pref  = False
                shr_pref = True

                src, shr_size = generate_tfmm_lines(arr, ndims=3, nvars=13,
                                    dtype=a.dtype, soasz=32, shr_max=92190,
                                    platform='cuda_S_tensor_lines_15')
                block = (192, 1, 1)
                grid = (ceil(nele*p*p/block[0]), 1, 1)
            if p == 3:
                l1_pref  = False
                shr_pref = True
                shr_max = 60416
                opargs = {'ld_opt': True, 'st_opt': False, 'intl_opt': False,
                        'pipe_opt': False, 'max_coag': 0, 'compute_size': 0,
                        'shr_op_order': 'grs', 'shr_bdc': True}
                
                src, shr_size = generate_tfmm_managed(arr, ndims=3, nvars=13, 
                                    soasz=32, dtype=a.dtype, opargs=opargs,
                                    platform='cuda_tfmm_managed', flux='hyperns',
                                    shared_max=shr_max, ufc_size=4*(2**15), 
                                    block_dim=128, source_term=True)
                block = (128, 1, 1)
                grid = (ceil(nele*p/block[0]), 1, 1)
            if p == 2:
                l1_pref  = False
                shr_pref = True
                shr_max = 43008
                opargs = {'ld_opt': True, 'st_opt': False, 'intl_opt': False,
                        'pipe_opt': False, 'max_coag': 0, 'compute_size': 0,
                        'shr_op_order': 'grs', 'shr_bdc': True}
                src, shr_size = generate_tfmm_managed(arr, ndims=3, nvars=13, 
                                    soasz=32, dtype=a.dtype, opargs=opargs,
                                    platform='cuda_tfmm_managed', flux='hyperns',
                                    shared_max=shr_max, ufc_size=4*(2**15), 
                                    block_dim=128, source_term=True)
                block = (128, 1, 1)
                grid = (ceil(nele*p/block[0]), 1, 1)

        print(f'{nele=}, {shr_size=}, {p=}, {block=}, {grid=}')

        f = open(f'tmul_S_p{p-1}.cu', 'w')
        f.write(src)
        f.close()

        # Build
        fun = self._build_kernel('gimmik_tfmm', src,
                                 [np.int32, np.intp]*2 + [np.int32])
        fun.set_cache_pref(prefer_l1=l1_pref, prefer_shared=shr_pref)
        fun.set_shared_size(dynm_shared=shr_size)

        class MulKernel(ComputeKernel):
            def run(self, queue):
                fun.exec_async(grid, block, queue.stream_comp, nele, b,
                               b.leaddim, out, out.leaddim,
                               dynm_shared=shr_size)

        return MulKernel()

