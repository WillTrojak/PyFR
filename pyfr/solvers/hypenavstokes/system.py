# -*- coding: utf-8 -*-

from pyfr.solvers.hypenavstokes.elements import HypeNavierStokesElements
from pyfr.solvers.hypenavstokes.inters import (HypeNavierStokesIntInters, HypeNavierStokesMPIInters,
                                         HypeNavierStokesBaseBCInters)
from pyfr.solvers.baseadvec import BaseAdvectionSystem


class HypeNavierStokesSystem(BaseAdvectionSystem):
    name = 'hyper-ns'

    elementscls = HypeNavierStokesElements
    intinterscls = HypeNavierStokesIntInters
    mpiinterscls = HypeNavierStokesMPIInters
    bbcinterscls = HypeNavierStokesBaseBCInters
