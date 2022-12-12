# -*- coding: utf-8 -*-

import numpy as np

from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionElements
from pyfr.solvers.mpeuler.elements import BaseMPFluidElements


class MPNavierStokesElements(BaseMPFluidElements, 
                             BaseAdvectionDiffusionElements):
    # Use the density field for shock sensing
    shockvar = 'rho'

    @property
    def _scratch_bufs(self):

        bufs = {'scal_fpts', 'vect_fpts', 'vect_upts', 'vect_upts_cpy',
                'scal_upts_cpy'}

        if 'flux' in self.antialias:
            bufs |= {'vect_qpts_cpy', 'scal_qpts', 'vect_qpts'}
        elif 'fraction' in self.antialias:
            bufs |= {'scal_qpts'}

        return bufs

    @staticmethod
    def grad_con_to_pri(cons, grad_cons, cfg):
        ns = cfg.getint('solver', 'species')

        arho, grad_arho = cons[:ns], grad_cons[:ns]
        rhouvw, grad_rhouvw = cons[ns:-ns], grad_cons[ns:-ns]
        rho, grad_rho = sum(arho), sum(grad_arho)

        E, grad_E = cons[-ns], grad_cons[-ns]

        alpha = np.vstack((cons[1-ns:], [1 - sum(cons[1-ns:])]))
        grad_a = np.vstack((grad_cons[1-ns:], [-sum(grad_cons[1-ns:])]))

        # Divide momentum components by ρ
        uvw = [rhov/rho for rhov in rhouvw]

        # Compute the specific energy
        gamma = [cfg.getfloat('constants', f'gamma{i}') for i in range(ns)]
        rhoe = (E - 0.5*rho*sum(v*v for v in uvw))

        # Velocity gradients: ∇u⃗ = 1/ρ·[∇(ρu⃗) - u⃗ ⊗ ∇ρ]
        grad_uvw = [(grad_rhov - v*grad_rho) / rho
                    for grad_rhov, v in zip(grad_rhouvw, uvw)]

        # Phase density gradients: grad_Rho = (grad_arho - rho*grad_a)/a
        with np.errstate(divide='ignore', invalid='ignore'):
            grad_Rho = np.where(alpha[:,None,...] != 0, 
                (grad_arho - rho*grad_a)/alpha[:,None,...], 0)

        grad_p = grad_E - 0.5*(np.einsum('ijk,iljk->ljk', uvw, grad_rhouvw) +
                               np.einsum('ijk,iljk->ljk', rhouvw, grad_uvw))
        inv_agm = 1/sum(alpha[i]/(gamma[i] - 1) for i in range(ns))
        agm_grad = sum(grad_a[i]/(gamma[i] - 1) for i in range(ns))
        agm_grad *= inv_agm
        grad_p += rhoe*agm_grad
        grad_p *= inv_agm

        return np.vstack((grad_Rho, grad_uvw, [grad_p], grad_a[:ns-1]))

    def set_backend(self, backend, nscalupts, nonce, linoff):
        super().set_backend(backend, nscalupts, nonce, linoff)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.mpnavstokes.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.mpnavstokes.kernels.tfluxlin')

        # Handle shock capturing and Sutherland's law
        shock_capturing = self.cfg.get('solver', 'shock-capturing')
        visc_corr = self.cfg.get('solver', 'viscosity-correction', 'none')
        if visc_corr not in {'sutherland', 'none'}:
            raise ValueError('Invalid viscosity-correction option')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nspec': self.nspec,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs,
            'shock_capturing': shock_capturing,
            'visc_corr': visc_corr
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat
        av = self.artvisc

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), f=s(self._vect_upts, c),
                artvisc=s(av, c), smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), f=s(self._vect_qpts, c),
                artvisc=s(av, c), smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), f=s(self._vect_upts, l),
                artvisc=s(av, l), verts=self.ploc_at('linspts', l),
                upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                artvisc=s(av, l), verts=self.ploc_at('linspts', l),
                upts=self.qpts
            )
