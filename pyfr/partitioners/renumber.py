from collections import defaultdict
import re

import numpy as np

from pyfr.inifile import Inifile
from pyfr.partitioners.base import BasePartitioner


class RenumberingPartitioner(BasePartitioner):
    name = 'renumber'
    has_multiple_constraints = False

    # Integer options
    int_opts, enum_opts, dflt_opts = {}

    def __init__(self):
        pass

    def partition(self, mesh, rnum, progress=NullProgressSequence):
        uuid = mesh['mesh_uuid']
        soln_regex = re.compile(r'rnum_([a-z]+)_p([0-9]+)$')

        def part_soln_fn(soln):
            newsoln = defaultdict(list)
            newsoln['mesh_uuid'] = uuid

            for k, v in rnum.array_info('rnum').items():
                etype, part = soln_regex.match(k).groups()

                # Inspect first to get size
                inew, pold, iold = v[0]
                nupts, nvars = soln[f'soln_{etype}_p{pold}'].shape()[:-1]
                psoln = np.zeros((nupts, nvars, len(v)), dtype=float)

                for inew, pold, iold in v:
                    psoln[..., inew] = soln[f'soln_{etype}_p{pold}'][..., iold]

                newsoln[f'soln_{etype}_p{part}'] = psoln

            # Copy previous config files
            newsoln |= rnum.array_info('config')

            return newsoln

        return part_soln_fn
