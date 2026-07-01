import numpy as np

from pyfr.backends.base.linalg import BaseGPULinalgKernels
from pyfr.backends.metal.provider import MetalKernel, MetalKernelProvider


class MetalLinalgKernels(BaseGPULinalgKernels, MetalKernelProvider):
    _sgsize = 32

    def _batched_inv(self, m, out, *, eidxs):
        n, neles = out.block_size, m.ioshape[2]
        cdtype = self._compute_dtype(m.dtype)
        ftsz = np.dtype(cdtype).itemsize

        inv_kern, nthreads = self._compile_inv(out, m.dtype)

        mem_alloc = self.backend.mem_alloc
        chunk = min(neles, 1024)
        scw = self._inv_sc(n, cdtype)
        anb, cnb = chunk*ftsz*n*n, chunk*ftsz*n*scw

        class BatchedInvKernel(MetalKernel):
            scale = 1.0

            def bind(self, *, scale=1.0):
                self.scale = scale

            def run(self, cbuf):
                sa, scbuf = mem_alloc(anb), mem_alloc(cnb)
                for eoff in range(0, neles, chunk):
                    b = min(chunk, neles - eoff)
                    inv_kern(cbuf, (b*nthreads, 1, 1), (nthreads, 1, 1),
                             (sa, 0), (scbuf, 0), out.data, m.data, eidxs.data,
                             m.leaddim, eoff, self.scale)

        return BatchedInvKernel(mats=[m, out, eidxs])

    def _batched_tiled_matvec(self, x, minv, y, tplargs):
        ixdtype = self.backend.ixdtype

        tpl = self.backend.lookup.get_template('batched_tiled_matvec')

        kern = self._build_kernel(
            'batched_tiled_matvec', tpl.render(**tplargs),
            [np.uintp, np.uintp, ixdtype, np.uintp, ixdtype]
        )

        nt = self._mv_nthreads
        grid, tgrp = (minv.nmats*nt, 1, 1), (nt, 1, 1)
        kargs = [minv.data, x.data, x.leaddim, y.data, y.leaddim]

        class ApplyTiledMatrixKernel(MetalKernel):
            def run(self, cbuf):
                kern(cbuf, grid, tgrp, *kargs)

        return ApplyTiledMatrixKernel(mats=[x, minv, y])
