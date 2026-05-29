from functools import cached_property

import numpy as np

from pyfr.shapes import BaseShape
from pyfr.util import subclass_where


def split_samples(samples, nvars):
    # Primitives are the first nvars rows; remaining rows form grad blocks
    pris = list(samples[:nvars])
    if samples.shape[0] > nvars:
        return pris, np.split(samples[nvars:], nvars)
    else:
        return pris, None


class VolumePostProcData:
    def __init__(self, cfg, pris, grad_pris=None):
        self.cfg = cfg
        self.pris = pris
        self.grad_pris = grad_pris
        self.fields = {}

    @property
    def nvars(self):
        return len(self.pris)

    @property
    def has_grads(self):
        return self.grad_pris is not None


class BoundaryPostProcData(VolumePostProcData):
    def __init__(self, cfg, pris, spts, etype, fidx, svpts, grad_pris=None):
        super().__init__(cfg, pris, grad_pris=grad_pris)
        self._spts = spts
        self._etype = etype
        self._fidx = fidx
        self._svpts = svpts

    @cached_property
    def _elementscls(self):
        from pyfr.solvers.base import BaseSystem
        sname = self.cfg.get('solver', 'system')
        return subclass_where(BaseSystem, name=sname).elementscls

    @cached_property
    def _shape(self):
        cls = subclass_where(BaseShape, name=self._etype)
        return cls(len(self._spts), self.cfg)

    @cached_property
    def _eles(self):
        return self._elementscls(type(self._shape), self._spts, self.cfg)

    @cached_property
    def pnorm(self):
        _, _, norm = self._shape.faces[self._fidx]
        norm_tiled = np.tile(norm, (len(self._svpts), 1))
        pn = self._eles.pnorm_at(self._svpts, norm_tiled)

        return pn.transpose(2, 0, 1)

    @cached_property
    def normals(self):
        return self.pnorm / np.linalg.norm(self.pnorm, axis=0)

    @cached_property
    def min_upt_wall_dist_approx(self):
        shape = self._shape
        _, proj, norm = shape.faces[self._fidx]
        upts = shape.upts

        norm /= np.linalg.norm(norm)
        face_pt = proj(*([0]*(shape.ndims - 1)))
        t = (face_pt - upts) @ norm
        upts_on_face = upts + t[:, None] * norm

        sbasis = shape.sbasis

        def interp(op):
            r = op @ self._spts.reshape(op.shape[1], -1)
            return r.reshape(op.shape[0], *self._spts.shape[1:])

        x_upt = interp(sbasis.nodal_basis_at(upts))
        x_face = interp(sbasis.nodal_basis_at(upts_on_face))

        dist = np.linalg.norm(x_upt - x_face, axis=2)

        return dist[t != 0].min(axis=0)
