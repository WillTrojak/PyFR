# -*- coding: utf-8 -*-

from itertools import chain

from pyfr.integrators.dual.phys.base import BaseDualIntegrator


class BaseDualStepper(BaseDualIntegrator):

    def _set_stage_n(self, s):
        pass

    def _finalize_stage(self, tcurr, s):
        pass


class BaseBDFStepper(BaseDualStepper):
    @property
    def _nstages(self):
        return 1

    @property
    def _stage_nregs(self):
        return 0

    @property
    def _stepper_nregs(self):
        return len(self._stepper_static_coeffs) - 1

    def _finalize_stage(self, tcurr, s):
        psnregs = self.pseudointegrator._pseudo_stepper_nregs

        # Rotate the source registers to the right by one
        self.pseudointegrator._regidx[
            psnregs:psnregs + self.pseudointegrator._stepper_nregs
        ] = (self.pseudointegrator._stepper_regidx[-1:]
             + self.pseudointegrator._stepper_regidx[:-1])

        # Copy the current soln into the first source register
        self.pseudointegrator._add(
            0, self.pseudointegrator._stepper_regidx[0],
            1, self.pseudointegrator._idxcurr
        )

    @property
    def _stepper_coeffs(self):
        return [1] + [sc/self._dt for sc in self._stepper_static_coeffs]


class DualBDF2Stepper(BaseBDFStepper):
    stepper_name = 'bdf2'

    @property
    def _stepper_order(self):
        return 2

    @property
    def _stepper_static_coeffs(self):
        return [-1.5, 2.0, -0.5]


class DualBDF3Stepper(BaseBDFStepper):
    stepper_name = 'bdf3'

    @property
    def _stepper_order(self):
        return 3

    @property
    def _stepper_static_coeffs(self):
        return [-11.0/6.0, 3.0, -1.5, 1.0/3.0]


class DualBackwardEulerStepper(BaseBDFStepper):
    stepper_name = 'backward-euler'

    @property
    def _stepper_order(self):
        return 1

    @property
    def _stepper_static_coeffs(self):
        return [-1.0, 1.0]


class BaseDIRKStepper(BaseDualStepper):

    @property
    def _stepper_nregs(self):
        return 1

    @property
    def _stage_nregs(self):
        return self._nstages

    def _set_stage_n(self, s):
        self._current_stage = s

    def _finalize_stage(self, tcurr, s):
        # call system.rhs with idxcurr and store it somewhere
        # store the time derivative of current stage for later
        self.pseudointegrator.system.rhs(
            tcurr, self.pseudointegrator._idxcurr,
            self.pseudointegrator._stage_regidx[s]
        )
        # finalize dirk
        if s == self._nstages - 1:
            #self.pseudointegrator._add(
            #    0, self.pseudointegrator._idxcurr,
            #    1, self.pseudointegrator._stepper_regidx[0],
            #    *chain(*zip([bred*self._dt for bred in self.b],
            #                self.pseudointegrator._stage_regidx))
            #)
            # if b[:] == a[_nstage][:]
            # last idxcurr is the new solution

            # Copy the current soln into the first source register
            self.pseudointegrator._add(
                0, self.pseudointegrator._stepper_regidx[0],
                1, self.pseudointegrator._idxcurr
            )

    @property
    def _get_current_stage_n(self):
        return self._current_stage

    @property
    def _stepper_coeffs(self):
        csn = self._get_current_stage_n
        return [self.a[csn][csn]] \
            + [c/self._dt for c in [-1.0, 1.0]] \
            + self.a[csn][:csn]


class ESDIRK3Stepper(BaseDIRKStepper):
    stepper_name = 'esdirk3'

    a = [[], [], [], []]
    a[0] = [0]
    a[1] = [1767732205903/4055673282236, 1767732205903/4055673282236]
    a[2] = [2746238789719/10658868560708, -640167445237/6845629431997,
            1767732205903/4055673282236]
    a[3] = [1471266399579/7840856788654, -4482444167858/7529755066697,
            11266239266428/11593286722821, 1767732205903/4055673282236]

    b = a[3]

    @property
    def _nstages(self):
        return 4


class ESDIRK4Stepper(BaseDIRKStepper):
    stepper_name = 'esdirk4'

    @property
    def _nstages(self):
        return 6


class DIRKTestStepper(BaseDIRKStepper):
    stepper_name = 'dirktest'

    a = [[1]]

    b = a[0]

    @property
    def _nstages(self):
        return 1
