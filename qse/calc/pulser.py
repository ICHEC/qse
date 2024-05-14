"""
This module Provides QSE calculator interface to Pulser. It is derived from
the ASE's CP2K calculator, though at the end it may look very different from
ASE calculator.
https://pulser.readthedocs.io/en/stable/
"""


import os
import os.path
from warnings import warn
from subprocess import Popen, PIPE
import numpy as np
import ase.io
from qse.calc.calculator import (Calculator, all_changes, Parameters, CalculatorSetupError)
# from ase.calculators.calculator import (Calculator, all_changes, Parameters, CalculatorSetupError)

import pulser
from pulser_simulation import Simulation, SimConfig, QutipEmulator

class Pulser(Calculator):
    """QSE-Calculator for pulser.

    Pulser is an open-source Python software package.
    It provides easy-to-use libraries for designing and
    simulating pulse sequences that act on programmable
    arrays of neutral qbits, a promising platform for
    quantum computation and simulation.
    Online documentation: https://pulser.readthedocs.io
    White paper: Quantum 6, 629 (2022)
    Source code repository (go here for the latest docs): https://github.com/pasqal-io/Pulser
    License: Apache 2.0
    – see [LICENSE](https://github.com/pasqal-io/Pulser/blob/master/LICENSE) for details

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

    implemented_properties = ['energy', 'state', 'fidality']
    default_parameters = dict(
        auto_write=False,
        folder='./',
        max_scf=50,
        print_level='LOW')

    def __init__(self, restart=None, qbits=None, pulse=None, emulator=None,
                 label='pulser-run', debug=False, **kwargs):
        """Construct pulser-calculator object."""
        self.pulse = pulser.Pulse.ConstantPulse(400, 5*np.pi, 0, 0) if pulse is None else pulse
        self.emulator = QutipEmulator if emulator is None else emulator
        self._debug = debug
        self.label = None
        self.parameters = None
        self.results = None
        self.qbits = None

        Calculator.__init__(self, restart=restart, label=label, qbits=qbits, **kwargs)

        if restart is not None:
            self.read(restart)

    def __del__(self):
        """deleting process. Empty"""
        pass

    def set(self, **kwargs):
        """Set parameters like set(key1=value1, key2=value2, ...)."""
        msg = '"%s" is not a known keyword for the CP2K calculator. ' \
              'To access all features of CP2K by means of an input ' \
              'template, consider using the "inp" keyword instead.'
        for key in kwargs:
            if key not in self.default_parameters:
                raise CalculatorSetupError(msg % key)

        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write(self, label):
        'Write qbits, parameters and calculated results into restart files.'
        #if self._debug:
        #    print("Writing restart to: ", label)
        #self.qbits.write(label + '_restart.traj')
        #self.parameters.write(label + '_params.ase')
        #from ase.io.jsonio import write_json
        #with open(label + '_results.json', 'w') as fd:
        #    write_json(fd, self.results)
        pass

    def read(self, label):
        'Read qbits, parameters and calculated results from restart files.'
        #self.qbits = ase.io.read(label + '_restart.traj')
        #self.parameters = Parameters.read(label + '_params.ase')
        #from ase.io.jsonio import read_json
        #with open(label + '_results.json') as fd:
        #    self.results = read_json(fd)
        pass

    def calculate(self, device=None, progress=True):
        """Do the calculation.
        # system_changes=all_changes -> check it's relevance.
        """
        coords = self.qbits.positions[:,:2]
        reg = pulser.Register.from_coordinates(coords, prefix="q")
        dev = pulser.devices.MockDevice if device == None else device
        # we need to have/add an attribute to calc for device
        sequence = pulser.Sequence(reg, dev)
        sequence.declare_channel("ch0", "rydberg_global")
        sequence.add(self.pulse, "ch0")
        sequence.measure("ground-rydberg")
        sim = self.emulator.from_sequence(sequence)
        # check whether all the emulator `classes` have .from_sequence method.
        self.results = sim.run(progress_bar=progress)
        self.sequence = sequence
        self.sim = sim
        self.qbits.state = self.results.get_final_state()
        if self.parameters.auto_write:
            self.write(self.label)
    #
