# -*- coding: utf-8 -*-

import csv

from pyfr.integrators.dual.phys.base import BaseDualIntegrator


class BaseDualController(BaseDualIntegrator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Solution filtering frequency
        self._fnsteps = self.cfg.getint('soln-filter', 'nsteps', '0')

        self.step_file = self.cfg.get('solver-time-integrator', 'step-file', None)
        if self.step_file != 'None':
            with open(self.step_file, newline='') as f:
                reader = csv.reader(f, delimiter=',')
                self.p_nsteps = list(reader)
        else:
            print('Not using step file')


    def _accept_step(self, idxcurr):
        self.tcurr += self._dt
        self.nacptsteps += 1
        self.nacptchain += 1

        # Filter
        if self._fnsteps and self.nacptsteps % self._fnsteps == 0:
            self.pseudointegrator.system.filt(idxcurr)

        # Invalidate the solution cache
        self._curr_soln = None

        # Invalidate the solution gradients cache
        self._curr_grad_soln = None

        # Fire off any event handlers
        self.completed_step_handlers(self)

        # Abort if plugins request it
        self._check_abort()

        # Clear the pseudo step info
        self.pseudointegrator.pseudostepinfo = []


class DualNoneController(BaseDualController):
    controller_name = 'none'

    def advance_to(self, t):
        if t < self.tcurr:
            raise ValueError('Advance time is in the past')

        if self.step_file != 'None':
            for step in iter(self.p_nsteps):
                self.pseudointegrator.pseudo_advance(self.tcurr, nstep=int(step[1]))
                self._accept_step(self.pseudointegrator._idxcurr)
        else:
            while self.tcurr < t:
                self.pseudointegrator.pseudo_advance(self.tcurr)
                self._accept_step(self.pseudointegrator._idxcurr)