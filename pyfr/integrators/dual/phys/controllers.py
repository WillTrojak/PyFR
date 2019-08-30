# -*- coding: utf-8 -*-

from itertools import chain

from pyfr.integrators.dual.phys.base import BaseDualIntegrator


class BaseDualController(BaseDualIntegrator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Solution filtering frequency
        self._fnsteps = self.cfg.getint('soln-filter', 'nsteps', '0')

    def _accept_step(self, idxcurr):
        self.tcurr += self._dt
        self.nacptsteps += 1
        self.nacptchain += 1

        # Filter
        if self._fnsteps and self.nacptsteps % self._fnsteps == 0:
            self.pseudoitegrator.system.filt(idxcurr)

        # Invalidate the solution cache
        self._curr_soln = None

        # Fire off any event handlers
        self.completed_step_handlers(self)

        # Clear the pseudo step info
        self.pseudointegrator.pseudostepinfo = []


class DualNoneController(BaseDualController):
    controller_name = 'none'

    def advance_to(self, t):
        if t < self.tcurr:
            raise ValueError('Advance time is in the past')

        while self.tcurr < t:
            for s in range(self._nstages):
                print('stage n =', s, 't = ', self.tcurr)
                final_stage = True if s == self._nstages - 1 else False
                if self._nstages > 1:
                    self._initialize_dirk() if s == 0 \
                                            else self._advance_dirk_stage()
                self.pseudointegrator.pseudo_advance(
                    self.tcurr, self._stepper_coeffs, final_stage
                )
                # call system.rhs with idxcurr and store it somewhere
                if self._nstages > 1:
                    # store the time derivative of current stage for later
                    self.pseudointegrator.system.rhs(
                        self.tcurr, self.pseudointegrator._idxcurr,
                        self.pseudointegrator._stage_regidx[s]
                    )
                    # finalize dirk
                    if s == self._nstages - 1:
                        self.pseudointegrator._add(
                            0, self.pseudointegrator._idxcurr,
                            1, self.pseudointegrator._stepper_regidx[-1],
                            *chain(*zip(self.b,
                                        self.pseudointegrator._stage_regidx))
                        )
                        # if b[:] == a[_nstage][:]
                        # last idxcurr is the new solution
            self._accept_step(self.pseudointegrator._idxcurr)
