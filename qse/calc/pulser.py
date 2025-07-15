"""
This module Provides QSE calculator interface to Pulser. It is derived from
the ASE's CP2K calculator, though at the end it may look very different from
ASE calculator.
https://pulser.readthedocs.io/en/stable/
"""

from time import time

from qse import Signal
from qse.calc.calculator import Calculator

try:
    import pulser
    import pulser.waveforms
    import qutip
    from pulser_simulation import QutipEmulator

    CALCULATOR_AVAILABLE = True
except ImportError:
    CALCULATOR_AVAILABLE = False


class Pulser(Calculator):
    r"""
    QSE-Calculator for pulser.

    Parameters
    ----------
    amplitude: qse.Signal, pulser.waveforms.Waveform
        The amplitude pulse.
    detuning: qse.Signal, pulser.waveforms.Waveform
        The detuning pulse.
    qbits: qse.Qbits
        The qbits object.
    device
        Defaults to pulser.devices.MockDevice.
    emulator
        Defaults to QutipEmulator.
    label: str
        Defaults to "pulser-run".
    wtimes: bool
        Defaults to True.

    Notes
    -----
    Pulser is an open-source Python software package.
    It provides easy-to-use libraries for designing and
    simulating pulse sequences that act on programmable
    arrays of neutral qbits, a promising platform for
    quantum computation and simulation.
    Online documentation: https://pulser.readthedocs.io
    White paper: Quantum 6, 629 (2022)
    Source code: https://github.com/pasqal-io/Pulser
    License: Apache 2.0
    - see [LICENSE](https://github.com/pasqal-io/Pulser/blob/master/LICENSE)
    for details.
    """

    implemented_properties = ["energy", "state", "fidality"]

    def __init__(
        self,
        amplitude,
        detuning,
        qbits,
        device=None,
        emulator=None,
        label="pulser-run",
        wtimes=True,
    ):
        installation_message = (
            "Pulser is not installed. To install, "
            "see https://pulser.readthedocs.io/en/stable/installation.html."
        )

        super().__init__(
            CALCULATOR_AVAILABLE, installation_message, label=label, qbits=qbits
        )
        self.device = pulser.devices.MockDevice if device is None else device
        self.emulator = QutipEmulator if emulator is None else emulator
        self.wtimes = wtimes
        self.results = None
        self.channel = "rydberg_global"

        self.amplitude = _format_pulse(amplitude)
        self.detuning = _format_pulse(detuning)

        self._sequence = None
        self._sim = None

        self.spins = None
        self.sij = None

        self.build_sequence()

    @property
    def coords(self):
        if self._qbits is None:
            return None
        return self._qbits.positions[:, :2]

    @property
    def register(self):
        return pulser.Register.from_coordinates(self.coords)

    @property
    def sequence(self):
        return self._sequence

    def build_sequence(self):
        self._sequence = pulser.Sequence(self.register, self.device)
        self._sequence.declare_channel("ch0", self.channel)
        self._sequence.add(
            pulser.Pulse(amplitude=self.amplitude, detuning=self.detuning, phase=0),
            "ch0",
        )
        self._sequence.measure("ground-rydberg")
        self._sim = self.emulator.from_sequence(self._sequence)

    @property
    def sim(self):
        return self._sim

    def __del__(self):
        """deleting process. Empty"""
        pass

    def write(self, label):
        """
        Write qbits, parameters and calculated results into restart files.
        Not yet implemented.

        Parameters
        ----------
        label : string
            used in filename
        """
        pass

    def read(self, label):
        """
        Read qbits, parameters and calculated results from restart files.
        Not yet implemented.

        Parameters
        ----------
        label : string
            used in filename
        """
        pass

    def calculate(self, progress=True):
        """
        Do the calculation.
        # system_changes=all_changes -> check it's relevance.
        """
        # we need to have/add an attribute to calc for device
        # check whether all the emulator `classes` have .from_sequence method.
        if self.wtimes:
            t1 = time()
        self.results = self.sim.run(progress_bar=progress)

        final_state = self.results.get_final_state()
        self.statevector = qutip.core.dimensions.to_tensor_rep(final_state).flatten()

        self.spins = self.get_spins()

        if self.wtimes:
            t2 = time()
            print(f"time in compute and simulation = {t2 - t1} s.")


def _format_pulse(pulse):
    if isinstance(pulse, pulser.waveforms.Waveform):
        return pulse
    elif isinstance(pulse, Signal):
        return pulser.waveforms.InterpolatedWaveform(
            duration=pulse.duration, values=pulse.values
        )
    else:
        raise Exception(
            "Pulses must be either `qse.Signal` or" " `pulser.waveforms.Waveform`."
        )
