# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.edeulerdual.elements import EDDEulerElements
from pyfr.solvers.edeuler.inters import (EDEulerIntInters, EDEulerMPIInters,
                                         EDEulerBaseBCInters)


class EDDEulerSystem(BaseAdvectionSystem):
    name = 'ed-dual-euler'

    elementscls = EDDEulerElements
    intinterscls = EDEulerIntInters
    mpiinterscls = EDEulerMPIInters
    bbcinterscls = EDEulerBaseBCInters
