from pkg_resources import resource_listdir, resource_string
import re


class BaseStoredNASAPoly:
    def __init__(self, species):
        self._rpaths = rpaths = resource_listdir(__name__, '.')

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
        self.t_form = float(metadata[1][-1])

        self.cp = []
        for l1, l2, l3 in zip(Tdata[0::3], Tdata[1::3], Tdata[2::3]):
            coeffs = [(i, float(x)) for i, x in enumerate(l2, start=-1)]
            coeffs.extend((i, float(x)) for i, x in enumerate(l3[:2], start=3))
            b = float(l3[-2], l3[-1])
            self.cp.append((float(l1[0]), float(l1[1]), coeffs, b))

def get_nasa_poly(species, T=None, order=None):
    pass