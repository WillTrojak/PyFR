from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem
from pyfr.solvers.mpnavstokes.elements import MPNavierStokesElements
from pyfr.solvers.mpnavstokes.inters import (MPNavierStokesBaseBCInters,
                                             MPNavierStokesIntInters,
                                             MPNavierStokesMPIInters)


class NavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'mp-navier-stokes'

    elementscls = MPNavierStokesElements
    intinterscls = MPNavierStokesIntInters
    mpiinterscls = MPNavierStokesMPIInters
    bbcinterscls = MPNavierStokesBaseBCInters
