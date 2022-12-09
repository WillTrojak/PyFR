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
        kern_path = 'pyfr.solvers.mpnavstokes.kernels.'
        self._be.pointwise.register(kern_path + 'tflux')
        self._be.pointwise.register(kern_path + 'tfluxlin')
        self._be.pointwise.register(kern_path + 'negdivconf_mp')
        self._be.pointwise.register(kern_path + 'negdivconf_mp_aa')
        self._be.pointwise.register(kern_path + 'fgrad_copy')
        self._be.pointwise.register(kern_path + 'frac_matdiv')

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
        kernel = self._be.kernel

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), grad=s(self._vect_upts_cpy, c),
                f=s(self._vect_upts, c), artvisc=s(av, c), 
                smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), grad=s(self._vect_qpts_cpy, c),
                f=s(self._vect_qpts, c), artvisc=s(av, c), 
                smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), grad=s(self._vect_upts_cpy, l),
                f=s(self._vect_upts, l), artvisc=s(av, l), 
                verts=self.ploc_at('linspts', l), upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), grad=s(self._vect_qpts_cpy, l),
                f=s(self._vect_qpts, l), artvisc=s(av, l), 
                verts=self.ploc_at('linspts', l), upts=self.qpts
            )

        # What the source term expressions (if any) are a function of
        plocsrc = self._ploc_in_src_exprs

        # Source term kernel arguments
        srctplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nspec': self.nspec,
            'srcex': self._src_exprs
        }

        # Transformed to physical divergence kernel + source term
        plocupts = self.ploc_at('upts') if plocsrc else None

        self.kernels['copy_soln'] = lambda uin: kernel(
            'copy', self._scal_upts_cpy, self.scal_upts[uin]
        )

        if 'fraction' in self.antialias:
            self._fraction_aa_kernels(nonce, srctplargs)
            self.kernels['negdivconf_mp'] = lambda fout: kernel(
                'negdivconf_mp_aa', tplargs=srctplargs, 
                dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout], 
                rcpdjac=self.rcpdjac_at('upts'), mat=self._scal_upts_ugrad,
                ploc=plocupts, u=self._scal_upts_cpy
            )
        else:
            self.kernels['negdivconf_mp'] = lambda fout: kernel(
                'negdivconf_mp', tplargs=srctplargs, 
                dims=[self.nupts, self.neles], tdivtconf=self.scal_upts[fout], 
                rcpdjac=self.rcpdjac_at('upts'), ploc=plocupts, 
                u=self._scal_upts_cpy, grad=self._vect_upts_cpy,
            )

        # Re-register kernels for gradient calcs
        if self.basis.order > 0:
            self.kernels['tgradpcoru_upts'] = lambda uin: kernel(
                'mul', self.opmat('M4 - M6*M0'), self.scal_upts[uin],
                out=self._vect_upts_cpy
            )

        self.kernels['tgradcoru_upts'] = lambda: kernel(
            'mul', self.opmat('M6'), self._vect_fpts.slice(0, self.nfpts),
            out=self._vect_upts_cpy, beta=float(self.basis.order > 0)
        )

        if c in r:
            self.kernels['gradcoru_upts_curved'] = lambda: kernel(
                'gradcoru', tplargs=tplargs, dims=[self.nupts, r[c]],
                gradu=s(self._vect_upts_cpy, c), 
                smats=self.curved_smat_at('upts'),
                rcpdjac=self.rcpdjac_at('upts', c)
            )

        if l in r:
            self.kernels['gradcoru_upts_linear'] = lambda: kernel(
                'gradcorulin', tplargs=tplargs, dims=[self.nupts, r[l]],
                gradu=s(self._vect_upts_cpy, l),
                upts=self.upts, verts=self.ploc_at('linspts', l)
            )

        def gradcoru_fpts():
            nupts, nfpts = self.nupts, self.nfpts
            vupts, vfpts = self._vect_upts_cpy, self._vect_fpts

            # Exploit the block-diagonal form of the operator
            muls = [kernel('mul', self.opmat('M0'),
                           vupts.slice(i*nupts, (i + 1)*nupts),
                           vfpts.slice(i*nfpts, (i + 1)*nfpts))
                    for i in range(self.ndims)]

            return self._be.unordered_meta_kernel(muls)

        self.kernels['gradcoru_fpts'] = gradcoru_fpts

        if 'flux' in self.antialias and self.basis.order > 0:
            def gradcoru_qpts():
                nupts, nqpts = self.nupts, self.nqpts
                vupts, vqpts = self._vect_upts_cpy, self._vect_qpts_cpy

                # Exploit the block-diagonal form of the operator
                muls = [self._be.kernel('mul', self.opmat('M7'),
                                        vupts.slice(i*nupts, (i + 1)*nupts),
                                        vqpts.slice(i*nqpts, (i + 1)*nqpts))
                        for i in range(self.ndims)]

                return self._be.unordered_meta_kernel(muls)

            self.kernels['gradcoru_qpts'] = gradcoru_qpts

    def _fraction_aa_kernels(self, nonce, tplargs):
        kernel = self._be.kernel
        # Allocate extra memory for fraction material derivative
        tags = {'align'}
        ext = nonce + 'vect_upts_fgrad'
        self._vect_upts_fgrad = self._be.matrix((self.ndims, self.nupts, 
                                                 self.nspec - 1, self.neles),
                                                 tags=tags, extent=ext)
        ext = nonce + 'vect_qpts_fgrad'
        self._vect_qpts_fgrad = self._be.matrix((self.ndims, self.nqpts, 
                                                 self.nspec - 1, self.neles),
                                                 tags=tags, extent=ext)
        ext = nonce + 'scal_upts_ugrad'
        self._scal_upts_ugrad = self._be.matrix((self.nupts, self.nspec - 1, 
                                                 self.neles),                          
                                                 tags=tags, extent=ext)
        ext = nonce + 'scal_qpts_ugrad'
        self._scal_qpts_ugrad = self._be.matrix((self.nqpts, self.nspec - 1,
                                                 self.neles),
                                                 tags=tags, extent=ext)

        self.kernels['fgrad_copy'] = lambda uin: kernel(
            'fgrad_copy', tplargs=tplargs, dims=[self.nupts, self.neles],
            grad=self._vect_upts_cpy, fgrad=self._vect_upts_fgrad
        )

        def fgradcoru_qpts():
            nupts, nqpts = self.nupts, self.nqpts
            vupts, vqpts = self._vect_upts_fgrad, self._vect_qpts_fgrad

            # Exploit the block-diagonal form of the operator
            muls = [self._be.kernel('mul', self.opmat('M7'),
                                    vupts.slice(i*nupts, (i + 1)*nupts),
                                    vqpts.slice(i*nqpts, (i + 1)*nqpts))
                    for i in range(self.ndims)]

            return self._be.unordered_meta_kernel(muls)
        self.kernels['fgradcoru_qpts'] = fgradcoru_qpts

        self.kernels['frac_matdiv'] = lambda uin: kernel(
            'frac_matdiv', tplargs=tplargs, dims=[self.nqpts, self.neles],
            fgrad=self._vect_qpts_fgrad, u=self._scal_qpts, 
            mat=self._scal_qpts_ugrad
        )

        self.kernels['matd_proj'] =  lambda uin: kernel(
            'mul', self.opmat('M8'), self._scal_qpts_ugrad,
            out=self._scal_upts_ugrad
        )
