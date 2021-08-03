# -*- coding: utf-8 -*-

from gimmik import generate_mm, generate_mm_split
import numpy as np

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

        order = self.backend.cfg.get('solver', 'order', 0)
        split = order + 1

        # Determine the grid/block
        block_s = (32, 1, 1)
        block = (32*split, 1, 1)
        grid = get_grid_for_block(block_s, b.ncol)

        # Generate
        src = generate_mm_split(arr, dtype=a.dtype, platform='cuda',
                                alpha=alpha, beta=beta, block_dim=block[0],
                                split=split)

        # Build
        fun = self._build_kernel('gimmik_mm', src,
                                 [np.int32, np.intp]*2 + [np.int32])
        fun.set_cache_pref(prefer_l1=True)

        class MulKernel(ComputeKernel):
            def run(self, queue):
                fun.exec_async(grid, block, queue.stream_comp, b.ncol, b,
                               b.leaddim, out, out.leaddim)

        return MulKernel()
