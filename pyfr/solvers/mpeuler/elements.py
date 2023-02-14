import numpy as np

from pyfr.solvers.baseadvec import BaseAdvectionElements
from pyfr.solvers.mpeuler.polynomials import BaseStoredNASAPoly, get_species


class BaseMPFluidElements:

    def __init__(self, basiscls, eles, cfg):
        super().__init__(basiscls, eles, cfg)

        self.nspec = self.cfg.getint('solver', 'nspecies')

    @classmethod
    def privarmap(cls, cfg, ndims):
        ns = cfg.getint('solver', 'nspecies')
        return (['u', 'v', 'w'][:ndims] + ['p'] +
                [f'a{i}rho{i}' for i in range(ns)])

    @classmethod
    def convarmap(cls, cfg, ndims):
        ns = cfg.getint('solver', 'nspecies')
        return (['rhou', 'rhov', 'rhow'][:ndims] + ['E'] +
                [f'C{i}' for i in range(ns)])

    @classmethod
    def dualcoeffs(cls, cfg, ndims):
        return cls.convarmap(cfg, ndims)

    @classmethod
    def visvarmap(cls, cfg, ndims):
        ns = cfg.getint('solver', 'nspecies')
        return ([(f'density_{i}', [f'a{i}rho{i}']) for i in range(ns)] +
                [('velocity', ['u', 'v', 'w'][:ndims]), ('pressure', ['p'])])

    @staticmethod
    def pri_to_con(pris, cfg):
        ns = cfg.getint('solver', 'nspecies')
        Tref = cfg.getfloat('species', 'Tref')

        # Read in species
        spec = []
        for key in cfg.items('species', prefix='spec-'):
            species = cfg.get('species', key)
            spec.append(get_species(species, Tref))

        arho, p = pris[-ns:], pris[-(ns + 1)]
        rho = np.sum(arho, axis=0)
        C = [ar / s.m_weight for ar, s in zip(arho, spec)]

        # Multiply velocity components by rho
        rhovs = [rho*v for v in pris[:-(ns + 1)]]

        # Compute the energy
        T = p / np.sum([c*s.R0 for c, s in zip(C, spec)], axis=0)
        rhoh = sum(r*s.R*np.polyval(s.Hr, T) for r, s in zip(arho, spec))
        E = np.array([rhoh - p + 0.5*rho*sum(v*v for v in pris[:-(ns + 1)])])

        return rhovs + [E] + C

    @staticmethod
    def con_to_pri(cons, cfg):
        ns = cfg.getint('solver', 'nspecies')
        Tref = cfg.getfloat('species', 'Tref')

        # Read in species
        spec = []
        for key in cfg.items('species', prefix='spec-'):
            species = cfg.get('species', key)
            spec.append(get_species(species, Tref))

        C, E = cons[-ns:], cons[-(ns + 1)]
        M = np.array([[s.m_weight] for s in spec])
        arho = [c * s.m_weight for c, s in zip(C, spec)]
        rho = np.sum(arho, axis=0)

        # Divide momentum components by rho
        vs = [rhov/rho for rhov in cons[:-(ns + 1)]]

        # Compute internal energy
        rhoe = E - 0.5*rho*sum(v*v for v in vs)

        # Initial linear guess for temperature
        Cv_ref = np.sum([s.Cv_ref*s.m_weight*c for c, s in zip(C, spec)], axis=0)
        T = rhoe / Cv_ref

        # Newton solve for temperature
        k_max = 5
        for _ in range(k_max):
            rhoe_dash = sum(s.m_weight*s.R*c*(np.polyval(s.Hr, T) - T)
                            for c, s in zip(C, spec))

            drhoedt = sum(s.m_weight*s.R*c*np.polyval(s.Cv, T)
                          for c, s in zip(C, spec))
            dT = (rhoe - rhoe_dash) / drhoedt
            T += dT
        p = (sum(s.m_weight*s.R*c*np.polyval(s.Hr, T) for c, s in zip(C, spec))
             - rhoe)

        return np.vstack((vs, [p], arho))

    @staticmethod
    def temperature(rhoe, spec, C, kmax=5):
        Cv_ref = sum(s.Cv_ref*s.m_weight*c for c, s in zip(C, spec))
        T = rhoe / Cv_ref

        for k in range(kmax):
            rhoe_dash = sum(s.m_weight*s.R*c*(np.polyval(s.Hr, T) - T)
                            for c, s in zip(C, spec))

            drhoedt = sum(s.m_weight*s.R*c*np.polyval(s.Cv, T)
                          for c, s in zip(C, spec))
            dT = (rhoe - rhoe_dash) / drhoedt
            T += dT
        return T

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
                'pyfr.solvers.euler.kernels.entropylocal'
            )
            self._be.pointwise.register(
                'pyfr.solvers.euler.kernels.entropyfilter'
            )

            # Template arguments
            eftplargs = {
                'ndims': self.ndims,
                'nupts': self.nupts,
                'nfpts': self.nfpts,
                'nvars': self.nvars,
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
            eftplargs['d_min'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'd-min', 1e-6)
            eftplargs['p_min'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'p-min', 1e-6)

            # Entropy tolerance
            eftplargs['e_tol'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'e-tol', 1e-6)

            # Hidden kernel parameters
            eftplargs['f_tol'] = self.cfg.getfloat('solver-entropy-filter',
                                                   'f-tol', 1e-4)
            eftplargs['ill_tol'] = self.cfg.getfloat('solver-entropy-filter',
                                                     'ill-tol', 1e-6)
            eftplargs['niters'] = self.cfg.getfloat('solver-entropy-filter',
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

        # Load in species
        Tref = self.cfg.getfloat('species', 'Tref')
        spec = []
        for k in self.cfg.items('species', prefix='spec-'):
            species = self.cfg.get('species', k)
            spec.append(get_species(species, Tref))

        s_dict = {'Rconst': BaseStoredNASAPoly.R0}
        for i, s in enumerate(spec):
            s_dict |= s.as_dict(suffix=i)

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nverts': len(self.basis.linspts),
            'nspec': self.nspec,
            'c': self.cfg.items_as('constants', float),
            's': s_dict,
            'jac_exprs': self.basis.jac_exprs
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), f=s(self._vect_upts, c),
                smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), f=s(self._vect_qpts, c),
                smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), f=s(self._vect_upts, l),
                verts=self.ploc_at('linspts', l), upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                verts=self.ploc_at('linspts', l), upts=self.qpts
            )
