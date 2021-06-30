# -*- coding: utf-8 -*-

from pyfr.solvers.achdnavstokes.elements import ACHDNavierStokesElements
from pyfr.solvers.achdnavstokes.inters import (ACHDNavierStokesIntInters, 
                                               ACHDNavierStokesMPIInters,
                                               ACHDNavierStokesBaseBCInters)
from pyfr.solvers.baseadvec import BaseAdvectionSystem


class ACHDNavierStokesSystem(BaseAdvectionSystem):
    name = 'ac-hd-navierstokes'

    elementscls = ACHDNavierStokesElements
    intinterscls = ACHDNavierStokesIntInters
    mpiinterscls = ACHDNavierStokesMPIInters
    bbcinterscls = ACHDNavierStokesBaseBCInters
