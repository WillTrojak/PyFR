# -*- coding: utf-8 -*-

import numpy as np
from collections import defaultdict

from pyfr.solvers.baseadvec import BaseAdvectionElements

class BaseMPFluidElements:

    def __init__(self, basiscls, eles, cfg):
        super().__init__(basiscls, eles, cfg)
        
        self.nspec = self.cfg.getint('solver', 'species')

    @classmethod
    def privarmap(cls, cfg, ndims):
        ns = cfg.getint('solver', 'species')
        m = defaultdict(lambda: None)
        m |= {2: ([f'a{i}rho{i}' for i in range(ns)] + ['u', 'v', 'p']),
              3: ([f'a{i}rho{i}' for i in range(ns)] + ['u', 'v', 'w', 'p']),
             }
        return m[ndims]

    @classmethod
    def convarmap(cls, cfg, ndims):
        ns = cfg.getint('solver', 'species')
        m = defaultdict(lambda: None)
        m |= {2: ([f'a{i}rho{i}' for i in range(ns)] + 
                  ['rhou', 'rhov', 'E']),
              3: ([f'a{i}rho{i}' for i in range(ns)] + 
                  ['rhou', 'rhov', 'rhow', 'E']),
             }
        return m[ndims]

    @classmethod
    def dualcoeffs(cls, cfg, ndims):
        return cls.convarmap(cfg, ndims)

    @classmethod
    def visvarmap(cls, cfg, ndims):
        ns = cfg.getint('solver', 'species')
        m = defaultdict(lambda: None)
        m |= {2: [(f'density_{i}', [f'a{i}rho{i}']) for i in range(ns)] +
                 [('velocity', ['u', 'v']),
                  ('pressure', ['p'])],
              3: [(f'density_{i}', [f'a{i}rho{i}']) for i in range(ns)] + 
                 [('velocity', ['u', 'v', 'w']),
                  ('pressure', ['p'])],
             }
        return m[ndims]

    @staticmethod
    def pri_to_con(pris, cfg):
        ns = cfg.getint('solver', 'species')
        p = pris[-1]
        arho = pris[:ns]
        rho = sum(arho)

        # Multiply velocity components by rho
        rhovs = [rho*c for c in pris[ns:-1]]

        # Compute the energy
        cp = sum(arho[i]*cfg.getfloat('constants', f'cp{i}') for i in range(ns))
        cv = sum(arho[i]*cfg.getfloat('constants', f'cv{i}') for i in range(ns))
        gamma = cp / cv

        E = p/(gamma - 1) + 0.5*rho*sum(c*c for c in pris[ns:-1])

        return np.vstack((arho, rhovs, [E]))

    @staticmethod
    def con_to_pri(cons, cfg):
        ns = cfg.getint('solver', 'species')
        E = cons[-1]
        arho = cons[:ns]
        rho = sum(arho)

        # Divide momentum components by rho
        vs = [rhov/rho for rhov in cons[ns:-1]]

        # Compute the pressure
        cp = sum(arho[i]*cfg.getfloat('constants', f'cp{i}') for i in range(ns))
        cv = sum(arho[i]*cfg.getfloat('constants', f'cv{i}') for i in range(ns))
        gamma = cp / cv
        p = (E - 0.5*rho*sum(v*v for v in vs))*(gamma - 1)

        return np.vstack((arho, vs, [p]))

    @staticmethod
    def validate_formulation(ctrl):
        shock_capturing = ctrl.cfg.get('solver', 'shock-capturing', 'none')
        if ctrl.formulation == 'dual' and shock_capturing == 'entropy-filter':
            raise ValueError('Entropy filtering not compatible with '
                             'dual time stepping.')

        ctrlvardt = ctrl.controller_has_variable_dt
        if ctrlvardt and shock_capturing == 'entropy-filter':
            raise ValueError('Entropy filtering not compatible with '
                             'adaptive time stepping.')

    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide shock-capturing at p = 0
        shock_capturing = self.cfg.get('solver', 'shock-capturing', 'none')

        # Modified entropy filtering method using specific physical
        # entropy (without operator splitting for Navier-Stokes)
        # doi:10.1016/j.jcp.2022.111501
        if shock_capturing == 'entropy-filter' and self.basis.order != 0:
            self._be.pointwise.register(
                'pyfr.solvers.mpeuler.kernels.entropylocal'
            )
            self._be.pointwise.register(
                'pyfr.solvers.mpeuler.kernels.entropyfilter'
            )

            # Template arguments
            eftplargs = {
                'ndims': self.ndims,
                'nupts': self.nupts,
                'nfpts': self.nfpts,
                'nvars': self.nvars,
                'nspec': self.nspec,
                'nfaces': self.nfaces,
                'c': self.cfg.items_as('constants', float),
                'order': self.basis.order
            }

            # Check to see if running collocated solution/flux points
            m0 = self.basis.m0
            mrowsum = np.max(np.abs(np.sum(m0, axis=1) - 1.0))
            if np.min(m0) < -1e-8 or mrowsum > 1e-8:
                raise ValueError('Entropy filter requires flux points to be a '
                                 'subset of solution points or a convex '
                                 'combination thereof.')

            # Minimum density/pressure constraints
            eftplargs['ad_min'] = self.cfg.getfloat('solver-entropy-filter',
                                                    'ad-min', 1e-6)
            eftplargs['p_min'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'p-min', 1e-6)

            # Minimuim and maximum fraction constraints
            eftplargs['a_min'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'a-min', 1e-6)
            eftplargs['a_max'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'a-max', 1 - 1e-6)

            # Entropy tolerance
            eftplargs['e_tol'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'e-tol', 1e-6)

            # Hidden kernel parameters
            eftplargs['f_tol']   = self.cfg.getfloat('solver-entropy-filter',
                                                     'f-tol', 1e-4)
            eftplargs['ill_tol'] = self.cfg.getfloat('solver-entropy-filter',
                                                     'ill-tol', 1e-6)
            eftplargs['niters']  = self.cfg.getfloat('solver-entropy-filter',
                                                     'niters', 20)

            # Precompute basis orders for filter
            ubdegs = self.basis.ubasis.degrees
            eftplargs['ubdegs'] = [int(max(dd)) for dd in ubdegs]
            eftplargs['order'] = self.basis.order

            # Compute local entropy bounds
            self.kernels['local_entropy'] = lambda uin: self._be.kernel(
                'entropylocal', tplargs=eftplargs, dims=[self.neles],
                u=self.scal_upts[uin], entmin_int=self.entmin_int
            )

            # Apply entropy filter
            self.kernels['entropy_filter'] = lambda uin: self._be.kernel(
                'entropyfilter', tplargs=eftplargs, dims=[self.neles],
                u=self.scal_upts[uin], entmin_int=self.entmin_int,
                vdm=self.vdm, invvdm=self.invvdm
            )


class MPEulerElements(BaseMPFluidElements, BaseAdvectionElements):
    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.mpeuler.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.mpeuler.kernels.tfluxlin')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nspec': self.nspec,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat
        kernel = self._be.kernel

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), f=s(self._vect_upts, c),
                smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), f=s(self._vect_qpts, c),
                smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), f=s(self._vect_upts, l),
                verts=self.ploc_at('linspts', l), upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                verts=self.ploc_at('linspts', l), upts=self.qpts
            )
