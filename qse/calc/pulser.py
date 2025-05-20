"""
This module Provides QSE calculator interface to Pulser. It is derived from
the ASE's CP2K calculator, though at the end it may look very different from
ASE calculator.
https://pulser.readthedocs.io/en/stable/
"""

from time import time

import numpy as np
import pulser
import pulser.waveforms
import qutip

# from pulser_simulation import Simulation, SimConfig, QutipEmulator
from pulser_simulation import QutipEmulator

import qse.calc.magnetic as magnetic
from qse import Signal
from qse.calc.calculator import Calculator, CalculatorSetupError
from qse.calc.messages import CalculatorSetupError

# from ase.calculators.calculator import (Calculator, all_changes, Parameters, CalculatorSetupError)


class Pulser(Calculator):
    """QSE-Calculator for pulser.

    Pulser is an open-source Python software package.
    It provides easy-to-use libraries for designing and
    simulating pulse sequences that act on programmable
    arrays of neutral qbits, a promising platform for
    quantum computation and simulation.
    Online documentation: https://pulser.readthedocs.io
    White paper: Quantum 6, 629 (2022)
    Source code: https://github.com/pasqal-io/Pulser
    License: Apache 2.0
    - see [LICENSE](https://github.com/pasqal-io/Pulser/blob/master/LICENSE) for details

    > TODO: if there is any config related to OMP_NUM_THREADS etc, place them here.

    Arguments:
    auto_write: bool
        Flag to enable the auto-write mode. If enabled the
        ``write()`` routine is called after every
        calculation, which mimics the behavior of the
        ``FileIOCalculator``. Default is ``False``.
    basis_set: str
        Name of the basis set to be use.
        The default is ``DZVP-MOLOPT-SR-GTH``.
    basis_set_file: str
        Filename of the basis set file.
        Default is ``BASIS_MOLOPT``.
        Set the environment variable $CP2K_DATA_DIR
        to enabled automatic file discovered.
    charge: float
        The total charge of the system.  Default is ``0``.
    cutoff: float
        The cutoff of the finest grid level.  Default is ``400 * Rydberg``.
    debug: bool
        Flag to enable debug mode. This will print all
        communication between the CP2K-shell and the
        CP2K-calculator. Default is ``False``.
    force_eval_method: str
        The method CP2K uses to evaluate energies and forces.
        The default is ``Quickstep``, which is CP2K's
        module for electronic structure methods like DFT.

    max_scf: int
        Maximum number of SCF iteration to be performed for
        one optimization. Default is ``50``.
    print_level: str
        PRINT_LEVEL of global output.
        Possible options are:
        DEBUG Everything is written out, useful for debugging purposes only
        HIGH Lots of output
        LOW Little output
        MEDIUM Quite some output
        SILENT Almost no output
        Default is 'LOW'
    """

    implemented_properties = ["energy", "state", "fidality"]
    default_parameters = dict(
        auto_write=False,
        folder="./",
        max_scf=50,
        print_level="LOW",
        label="q",
        qbits=None,
    )

    def __init__(
        self,
        qbits=None,
        amplitude=None,
        detuning=None,
        device=None,
        emulator=None,
        label="pulser-run",
        wtimes=True,
    ):
        """Construct pulser-calculator object.
        we need qubits, amplitude, detuning, device and emulator.
        """
        Calculator.__init__(self, label=label, qbits=qbits)
        self.device = pulser.devices.MockDevice if device == None else device
        self.emulator = QutipEmulator if emulator is None else emulator
        self.label = label
        self.wtimes = wtimes
        self.parameters = None
        self.results = None
        self.channel = "rydberg_global"
        if amplitude is not None:
            if isinstance(amplitude, (Signal, pulser.waveforms.Waveform)):
                self.amplitude = amplitude
            else:
                self.amplitude = Signal(amplitude)
        if detuning is not None:
            if isinstance(detuning, (Signal, pulser.waveforms.Waveform)):
                self.detuning = detuning
            else:
                self.detuning = Signal(detuning)
        #
        # assert self.amplitude.duration == self.detuning.duration
        self.duration = self.amplitude.duration
        wa = isinstance(self.amplitude, pulser.waveforms.Waveform)
        wd = isinstance(self.detuning, pulser.waveforms.Waveform)
        if wa:
            amp = self.amplitude
        else:
            amp = pulser.waveforms.InterpolatedWaveform(
                duration=self.duration, values=self.amplitude.values
            )
        if wd:
            det = self.detuning
        else:
            det = pulser.waveforms.InterpolatedWaveform(
                duration=self.duration, values=self.detuning.values
            )

        self.pulse = pulser.Pulse(amplitude=amp, detuning=det, phase=0)

        # self.pulse = pulser.Pulse(
        #    amplitude=pulser.waveforms.InterpolatedWaveform(
        #        duration=self.duration,
        #        values=self.amplitude.values),
        #        duration=self.duration,
        #    detuning=pulser.waveforms.InterpolatedWaveform(
        #        values=self.detuning.values), phase=0
        # )
        # pulser part which defines Hamiltonian parameters done. #
        self._register, self._sequence, self._sim = None, None, None
        self.qbits = qbits if qbits is not None else None

    # end of init #

    # defines properties #
    @property
    def qbits(self):
        return self._qbits

    @Calculator.qbits.setter
    def qbits(self, qbits, prefix="q"):
        if qbits is None:
            self._qbits, self._coords, self._register, self._sequence = (
                None,
                None,
                None,
                None,
            )
        else:
            self._qbits = qbits
            self._coords = qbits.positions[:, :2]
            self._register = pulser.Register.from_coordinates(
                self.coords, prefix=prefix
            )
            self.sequence = pulser.Sequence

    #

    @property
    def coords(self):
        return self._coords

    @property
    def register(self):
        return self._register

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        self._sequence = sequence(self.register, self.device)
        self._sequence.declare_channel("ch0", self.channel)
        self._sequence.add(self.pulse, "ch0")
        self._sequence.measure("ground-rydberg")
        self._sim = self.emulator.from_sequence(self._sequence)

    @property
    def sim(self):
        return self._sim

    def build_sequence(self):
        pass

    def __del__(self):
        """deleting process. Empty"""
        pass

    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2, ...)."""
        msg = (
            '"%s" is not a known keyword for the Pulser calculator. '
            "To access all features of Pulser by means of an input "
            'template, consider using the "inp" keyword instead.'
        )
        for key in kwargs:
            if key not in self.default_parameters:
                raise CalculatorSetupError(msg % key)

        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write(self, label):
        "Write qbits, parameters and calculated results into restart files."
        # if self._debug:
        #    print("Writing restart to: ", label)
        # self.qbits.write(label + '_restart.traj')
        # self.parameters.write(label + '_params.ase')
        # from ase.io.jsonio import write_json
        # with open(label + '_results.json', 'w') as fd:
        #    write_json(fd, self.results)
        pass

    def read(self, label):
        "Read qbits, parameters and calculated results from restart files."
        # self.qbits = ase.io.read(label + '_restart.traj')
        # self.parameters = Parameters.read(label + '_params.ase')
        # from ase.io.jsonio import read_json
        # with open(label + '_results.json') as fd:
        #    self.results = read_json(fd)
        pass

    def calculate(self, progress=True):
        """Do the calculation.
        # system_changes=all_changes -> check it's relevance.
        """
        # we need to have/add an attribute to calc for device
        # check whether all the emulator `classes` have .from_sequence method.
        if self.wtimes:
            t1 = time()
        self.results = self.sim.run(progress_bar=progress)
        # self.qbits.states = self.results.get_final_state()
        final_state = self.results.get_final_state()
        self.statevector = qutip.core.dimensions.to_tensor_rep(final_state).flatten()
        # if self.parameters.auto_write:
        #    self.write(self.label)
        self.spins = self.get_spins()
        self.sij = self.get_sij()
        # self.struc_fac = self.structure_factor_from_sij()
        if self.wtimes:
            t2 = time()
            print(f"time in compute and simulation = {t2 - t1} s.")

    #

    def get_spins(self):
        if self.results is None:
            self.calculate()
        #
        nqbits = len(self.qbits)
        ibasis = magnetic.get_basis(2**nqbits, nqbits)
        si = magnetic.get_spins(self.statevector, ibasis, nqbits)
        return si

    def get_sij(self):
        if self.results is None:
            self.calculate()
        #
        nqbits = len(self.qbits)
        ibasis = magnetic.get_basis(2**nqbits, nqbits)
        sij = magnetic.get_sisj(self.final_state, ibasis, nqbits)
        return sij

    def structure_factor_from_sij(self, L1, L2, L3):
        struc_fac = magnetic.structure_factor_from_sij(L1, L2, L3, self.qbits, self.sij)
        return struc_fac
