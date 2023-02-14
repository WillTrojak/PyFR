from functools import cache
from math import comb, cos, log, pi
from pkg_resources import resource_listdir, resource_string

import numpy as np


def bernstein(x, i, n):
    return comb(n, i)*pow(x, i)*pow(1 - x, n - i)


def cheblobatto_points(n, a=-1, b=1):
    return 0.5*np.array([(a - b)*cos(pi*j/(n - 1)) + (b + a)
                         for j in range(n)])


class BaseStoredNASAPoly:
    R0 = 8314.46261815324
    T_ref_nasa = 298.15

    def __init__(self, species, Trange, deg=5):
        self._rpaths = rpaths = resource_listdir(__name__, '.')

        self.species = species
        self.Trange = Trange

        for path in rpaths:
            if f'{species}.txt' == path:
                self.path = path
                break
        else:
            raise ValueError(f"Rule not found for species '{species}'")

        raw = resource_string(__name__, f'./{self.path}').decode()
        self.data = [x.split() for x in raw.split('\n')]
        self.metadata = metadata = self.data[:2]
        self.Tdata = Tdata = self.data[2:]

        self.n_T_range = int(metadata[1][0])
        self.m_weight = float(metadata[1][-2])
        self.R = self.R0 / self.m_weight
        self.t_form = float(metadata[1][-1])

        self._cp = []
        for l1, l2, l3 in zip(Tdata[0::3], Tdata[1::3], Tdata[2::3]):
            coeff = [(i, float(x)) for i, x in enumerate(l2, start=-2)]
            coeff.extend((i, float(x)) for i, x in enumerate(l3[:2], start=3))
            T0, T1 = float(l1[0]), float(l1[1])
            b1, b2 = float(l3[-2]), float(l3[-1])

            h_offset = 1000*float(l1[-1])/self.m_weight
            if T0 <= self.T_ref_nasa <= T1:
                h_ref = self.R*(self._base_poly_int(coeff, self.T_ref_nasa) + b1)
                h_end = (self.R*(self._base_poly_int(coeff, T1) + b1)
                         + h_offset - h_ref)
                H_std = h_offset - h_ref
            else:
                h_0 = self.R*(self._base_poly_int(coeff, T0) + b1)
                h_1 = self.R*(self._base_poly_int(coeff, T1) + b1)
                H_std = h_end - h_0
                h_end = H_std + h_1

            self._cp.append({'T0': T0, 'T1': T1, 'C': coeff, 'b1': b1, 'b2': b2,
                             'H_std': H_std, 'h_offset': h_offset})

        coeff = [cp for cp in self._cp if cp['T0'] <= Trange <= cp['T1']][0]

        # Reinterpolate NASA coefficent to polynomial
        self.reinterpolate(coeff, deg, ptype='poly')

        # Reference cp, cv and gamma
        self.Cp_ref = self.R*np.polyval(self.Cp, self.T_ref_nasa)
        self.Cv_ref = self.Cp_ref - self.R
        self.gamma_ref = self.Cp_ref / self.Cv_ref

        # Enthalpy polynomial
        self.Hr = np.polyint(self.Cp, m=1)
        self.Hr[-1] = (coeff['h_offset'] / self.R
                       - np.polyval(self.Hr, self.T_ref_nasa))

        # Cv / R = (Cp / R) - 1
        self.Cv = np.hstack((self.Cp[:-1], self.Cp[-1] - 1))

    def as_dict(self, suffix=None):
        s_str = str(suffix) if suffix is not None else ''
        params = [('Cp', self.Cp), ('Cv', self.Cv), ('Hr', self.Hr),
                  ('M', self.m_weight), ('R', self.R), ('Trange', self.Trange),
                  ('sepcies', self.species),
                  ('Cv_ref', self.Cv_ref), ('Cp_ref', self.Cv_ref),
                  ('gamma_ref', self.gamma_ref)]
        return {f"{k}{s_str}": v for k, v in params}

    @staticmethod
    def _base_poly_eval(p, x):
        z = 0*x
        for o, c in p:
            z += c*pow(x, o)
        return z

    @staticmethod
    def _base_poly_int(p, x):
        z = np.zeros_like(x)
        for o, c in p:
            if o == -1:
                z += c*log(x)
            else:
                z += c*pow(x, o + 1)/(o + 1)
        return z

    def reinterpolate(self, coeff, deg, ptype='bern'):
        self.Tmin, self.Tmax = coeff['T0'], coeff['T1']

        match ptype.lower():
            case 'bern':
                Ts_ref = cheblobatto_points(deg + 1, 0, 1)
                Ts = cheblobatto_points(deg + 1, coeff['T0'], coeff['T1'])
                Tb = np.linspace(coeff['T0'], coeff['T1'], deg + 1)

                Cp_b = self._base_poly_eval(coeff['C'], Tb)
                Cp_s = 0*Cp_b
                for i, b in enumerate(Cp_b):
                    Cp_s += b*bernstein(Ts_ref, i, deg)
            case 'poly':
                # Lobatto to ensure continuity between ranges
                Ts = cheblobatto_points(deg + 1, coeff['T0'], coeff['T1'])
                Cp_s = self._base_poly_eval(coeff['C'], Ts)
            case _:
                raise KeyError('NASA poly reinterpolation type not supported')

        self.Cp = np.polyfit(Ts, Cp_s, deg)

@cache
def get_species(species, Trange):
    return BaseStoredNASAPoly(species, Trange)
