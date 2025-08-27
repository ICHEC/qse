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
    device: pulser.devices.Device
        The device.
        Defaults to pulser.devices.MockDevice.
    emulator
        The emulator.
        Defaults to QutipEmulator.
    label: str
        The label.
        Defaults to "pulser-run".
    wtimes: bool
        Whether to print the times.
        Defaults to True.

    Examples
    --------
    A simple example of using the Pulser calculator:

    >>> qbits = qse.lattices.chain(4.0, 4)
    >>> duration = 400
    >>> pulser_calc = qse.calc.Pulser(
    ...     amplitude=qse.Signal(np.ones(6) * 1.01, duration),
    ...     detuning=qse.Signal(np.ones(6) * 0.12, duration),
    ...     qbits=qbits,
    ... )
    >>> pulser_calc.build_sequence()
    >>> pulser_calc.calculate()
    >>> pulser_calc.get_spins()

    Notes
    -----
    Pulser will only use the x-y coordinates of the Qbits object.

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
        qbits=None,
        amplitude=None,
        detuning=None,
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
            qbits=qbits, label=label, is_calculator_available=CALCULATOR_AVAILABLE, installation_message=installation_message, 
        )
        self.device = pulser.devices.MockDevice if device is None else device
        self.emulator = QutipEmulator if emulator is None else emulator
        self.wtimes = wtimes
        self.results = None
        self.channel = "rydberg_global"

        self.amplitude = amplitude
        self.detuning = detuning

        self._sequence = None
        self._sim = None
        self.statevector = None

    @property
    def amplitude(self):
        return self._amplitude

    @amplitude.setter
    def amplitude(self, amplitude):
        self._amplitude = _format_pulse(amplitude)

    @property
    def detuning(self):
        return self._detuning

    @detuning.setter
    def detuning(self, detuning):
        self._detuning = _format_pulse(detuning)

    @property
    def coords(self):
        if self._qbits is None:
            return None
        return self._qbits.positions[:, :2]

    @property
    def register(self):
        return pulser.Register.from_coordinates(self.coords, prefix="q")

    @property
    def sequence(self):
        return self._sequence

    def build_sequence(self):
        """
        Build the sequence of operations involving the qubit coordinates,
        amplitude pulse and detuning pulse.
        """
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
    if pulse is None or isinstance(pulse, pulser.waveforms.Waveform):
        return pulse
    if isinstance(pulse, Signal):
        return pulser.waveforms.InterpolatedWaveform(
            duration=pulse.duration, values=pulse.values
        )
    raise Exception(
        "Pulses must be either `qse.Signal` or `pulser.waveforms.Waveform`."
    )
