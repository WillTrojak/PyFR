# -*- coding: utf-8 -*-

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
                self._set_currstg(s)

                # Handle the first explicit step in ESDIRK
                if True: #not self._stepper_coeffs[0] == 0:
                    self.pseudointegrator.pseudo_advance(
                        self.tcurr, self._stepper_coeffs, self._currstg
                    )

            self._finalize_stage(self.tcurr)

            self._accept_step(self.pseudointegrator._idxcurr)
