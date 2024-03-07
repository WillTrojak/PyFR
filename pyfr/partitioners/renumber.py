from collections import defaultdict
import re

import numpy as np

from pyfr.inifile import Inifile
from pyfr.partitioners.base import BasePartitioner
from pyfr.progress import NullProgressSequence


class RenumberingPartitioner(BasePartitioner):
    name = 'renumber'
    has_multiple_constraints = True

    # Integer options
    int_opts, enum_opts, dflt_opts = {}, {}, {}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def partition(self, mesh, rnum, progress=NullProgressSequence):
        uuid = mesh['mesh_uuid']
        soln_regex = re.compile(r'rnum_([a-z]+)_p([0-9]+)$')

        parts = [n for n in rnum if n.startswith('rnum')]
        print(parts)

        def part_soln_fn(soln):
            newsoln = defaultdict(list)
            newsoln['mesh_uuid'] = uuid

            newsoln |= {n: soln[n] for n in soln if n.startswith('config')}

            for key in parts:
                etype, part = soln_regex.match(key).groups()

                idx = rnum[key]

                # Inspect first to get size
                inew, pold, iold = idx[0]

                nupts, nvars = soln[f'soln_{etype}_p{pold}'].shape[:-1]
                psoln = np.zeros((nupts, nvars, len(idx)), dtype=float)

                for inew, pold, iold in idx:
                    psoln[..., inew] = soln[f'soln_{etype}_p{pold}'][..., iold]

                newsoln[f'soln_{etype}_p{part}'] = psoln

            return newsoln

        return part_soln_fn
