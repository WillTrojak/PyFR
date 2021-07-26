# -*- coding: utf-8 -*-

from pyfr.solvers.edeuler.elements import EDEulerElements
from pyfr.solvers.edeuler.inters import (EDEulerIntInters, EDEulerMPIInters,
                                         EDEulerBaseBCInters)
from pyfr.solvers.baseadvec import BaseAdvectionSystem


class EDEulerSystem(BaseAdvectionSystem):
    name = 'ed-euler'

    elementscls = EDEulerElements
    intinterscls = EDEulerIntInters
    mpiinterscls = EDEulerMPIInters
    bbcinterscls = EDEulerBaseBCInters
