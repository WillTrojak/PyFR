from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.solvers.mpeuler.elements import MPEulerElements
from pyfr.solvers.mpeuler.inters import (MPEulerIntInters, MPEulerMPIInters,
                                         MPEulerBaseBCInters)


class MPEulerSystem(BaseAdvectionSystem):
    splitsystem = False
    name = 'mp-euler'

    elementscls = MPEulerElements
    intinterscls = MPEulerIntInters
    mpiinterscls = MPEulerMPIInters
    bbcinterscls = MPEulerBaseBCInters
