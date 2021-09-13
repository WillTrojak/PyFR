# -*- coding: utf-8 -*-

from pyfr.solvers.ednavstokes.inters import (EDNavierStokesBaseBCInters,
                                             EDNavierStokesIntInters,
                                             EDNavierStokesMPIInters)
from pyfr.solvers.ednavstokesdual.elements import EDDNavierStokesElements
from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem


class EDDNavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'ed-dual-navier-stokes'

    elementscls = EDDNavierStokesElements
    intinterscls = EDNavierStokesIntInters
    mpiinterscls = EDNavierStokesMPIInters
    bbcinterscls = EDNavierStokesBaseBCInters
