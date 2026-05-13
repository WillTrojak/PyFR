import os
from concurrent.futures import ThreadPoolExecutor
from functools import cached_property

import numpy as np

from pyfr.backends.base.linalg import (apply_tiled_tplargs,
                                       batched_inv_ok_dtypes)
from pyfr.backends.openmp.provider import OpenMPKernel, OpenMPKernelProvider
from pyfr.nputil import BLASThreadCtrl, bf16_to_f32, f32_to_bf16, is_bf16
from pyfr.util import ndrange


class OpenMPLinalgKernels(OpenMPKernelProvider):
    def __init__(self, backend):
        super().__init__(backend)

        # Capture the CPU allocation before OpenMP pins the master thread
        if hasattr(os, 'sched_getaffinity'):
            self.cpus = os.sched_getaffinity(0)
        else:
            self.cpus = None

    def _omp_num_threads(self):
        if raw := os.environ.get('OMP_NUM_THREADS'):
            return int(raw.split(',')[0])
        elif self.cpus:
            return len(self.cpus)
        else:
            return os.cpu_count()

    @cached_property
    def _inv_executor(self):
        cpus = self.cpus

        # Reset workers off the master's pinned core onto the full allocation
        init = (lambda: os.sched_setaffinity(0, cpus)) if cpus else None

        return ThreadPoolExecutor(max_workers=self._omp_num_threads(),
                                  thread_name_prefix='pyfr-omp-inv',
                                  initializer=init)

    def batched_inv_tiled(self, m, out, *, eidxs):
        if (not batched_inv_ok_dtypes(m.dtype, out.dtype) or
            out.block_size != m.nrow):
            raise ValueError('Invalid tiled output matrix')

        nrow, csubsz, nb = m.nrow, m.backend.csubsz, m.datashape[0]
        nworkers = min(self._omp_num_threads(), nb)
        executor = self._inv_executor if nworkers > 1 else None
        tile_ij = list(ndrange(out.ntiles, out.ntiles))
        ts = out.tsize
        ibf16, obf16 = is_bf16(m.dtype), is_bf16(out.dtype)
        eye = np.eye(nrow, dtype=np.float32 if ibf16 else m.dtype)

        earr = eidxs.get().ravel()
        blocks, rem = divmod(earr, out.csubsz)
        soas, inners = divmod(rem, out.soasz)

        class BatchedInvTiledKernel(OpenMPKernel):
            _scale = 1.0

            def bind(self, *, scale=1.0):
                self._scale = scale

            def run(self):
                data = m.data.reshape(m.datashape)

                def process(b0, b1):
                    start = b0*csubsz
                    end = min(b1*csubsz, len(earr))

                    jac = data[b0:b1].transpose(0, 2, 4, 1, 3)
                    jac = jac.reshape(-1, nrow, nrow)[:end - start]
                    if ibf16:
                        jac = bf16_to_f32(jac)
                    inv = np.linalg.inv(eye - self._scale*jac)

                    sl = slice(start, end)
                    b, s, n = blocks[sl], soas[sl], inners[sl]

                    for i, j in tile_ij:
                        r0, r1 = i*ts, min((i + 1)*ts, nrow)
                        c0, c1 = j*ts, min((j + 1)*ts, nrow)
                        blk = inv[:, r0:r1, c0:c1]
                        if obf16:
                            blk = f32_to_bf16(blk)
                        out.data[b, i, j, :r1 - r0, s, :c1 - c0, n] = blk

                if nworkers > 1:
                    with BLASThreadCtrl.serial():
                        # Submit the inversion work to the thread pool
                        fs = [executor.submit(process, i*nb//nworkers,
                                              (i + 1)*nb//nworkers)
                              for i in range(nworkers)]

                        # Wait for completion
                        for f in fs:
                            f.result()
                else:
                    process(0, nb)

        return BatchedInvTiledKernel(mats=[m, out, eidxs])

    def batched_tiled_matvec(self, x, minv, y, *, nupts, nvars,
                             in_scale=(), out_scale=()):
        tplargs = apply_tiled_tplargs(minv, nupts=nupts, nvars=nvars,
                                      in_scale=in_scale, out_scale=out_scale)

        tpl = self.backend.lookup.get_template('batched_tiled_matvec')
        src = tpl.render(**tplargs)

        ixdtype = self.backend.ixdtype
        argt = [ixdtype, np.uintp, np.uintp, np.uintp]
        kern = self._build_kernel('batched_tiled_matvec', src, argt)
        kern.set_args(minv.nmats, minv, x, y)

        return OpenMPKernel(mats=[x, minv, y], kernel=kern)
