# -*- coding: utf-8 -*-

from pyfr.solvers.ednavstokes.elements import EDNavierStokesElements
from pyfr.solvers.ednavstokes.inters import (EDNavierStokesBaseBCInters,
                                             EDNavierStokesIntInters,
                                             EDNavierStokesMPIInters)
from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem


class EDNavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'ed-navier-stokes'

    elementscls = EDNavierStokesElements
    intinterscls = EDNavierStokesIntInters
    mpiinterscls = EDNavierStokesMPIInters
    bbcinterscls = EDNavierStokesBaseBCInters
