"""
This module Provides QSE calculator interface to Pulser. It is derived from
the ASE's CP2K calculator, though at the end it may look very different from
ASE calculator.

It defines the `Pulser backend <https://pulser.readthedocs.io/en/stable/>`_
for analog computation.
"""

from time import time

from qse import Signal, Signals
from qse.calc.calculator import Calculator

try:
    import pulser
    import pulser.waveforms
    from pulser_simulation import QutipEmulator

    CALCULATOR_AVAILABLE = True
except ImportError:
    CALCULATOR_AVAILABLE = False

from functools import partial
from typing import Any, Dict, Generator

import numpy as np

from qse.calc.results import SimResult



# Function to extract metadata from Pulser's simulation
def extract_safe_metadata(obj, backend_name: str) -> dict:
    """
    Dynamically sweep an object for primitive attributes to build metadata.
    Ignore callables, complex objects, and dunder methods.
    TODO: Need to investigate some _ functions that give useful info.
    """
    metadata = {"backend": backend_name}
    
    for attr in dir(obj):
        # Optionally skip dunders, but Pulser keeps good stuff in _basis_name
        if attr.startswith('__'): 
            continue
            
        try:
            val = getattr(obj, attr)
            # Only extract basic serializable types
            if isinstance(val, (int, float, str, bool, tuple)):
                # Clean up the key name (e.g., '_basis_name' -> 'basis_name')
                clean_key = attr.lstrip('_')
                metadata[clean_key] = val
        except Exception:
            pass # Ignore properties that raise errors on access
    return metadata

# ==========================================
# GLOBAL EXTRACTORS (Pure Functions)
# ==========================================
# These functions know nothing about the Calculator.
# They only know how to translate a Pulser result into QSE format.


def extract_pulser_statevector(raw_results: Any) -> np.ndarray:
    """Extracts the final dense complex array."""
    qobj = raw_results.get_final_state()
    state_array = qobj.full().flatten()
    # this is equivalent to earlier flipping based on
    # the channel name == "rydberg_global"
    if getattr(raw_results, "_basis_name", "") == "ground-rydberg":
        state_array = state_array[::-1]
    return state_array


def extract_pulser_counts(raw_results: Any, shots: int) -> Dict[str, int]:
    """Samples the final state hardware style."""
    return dict(raw_results.sample_final_state(N_samples=shots))


def extract_pulser_expectation(raw_results: Any, observable: Any) -> float:
    """Delegates expectation value calculations to QuTiP."""
    exp_array = raw_results.expect([observable])[0]
    return float(exp_array[-1])


def extract_pulser_states(raw_results: Any) -> Generator[np.ndarray, None, None]:
    """Yields the statevector at every intermediate time step."""
    for qobj in raw_results.states:
        yield qobj.full().flatten()


# The calculator
class Pulser(Calculator):
    r"""
    QSE-Calculator for pulser.

    Parameters
    ----------
    qbits: qse.Qbits
        The qbits object.
    amplitude: qse.Signal, qse.Signals, pulser.waveforms.Waveform
        The amplitude pulse.
    detuning: qse.Signal, qse.Signals, pulser.waveforms.Waveform
        The detuning pulse.
    channel : str
        Which channel to use. For example "rydberg_global" for Rydberg or
        "mw_global" for microwave.
        Defaults to "rydberg_global".
    magnetic_field : np.ndarray | list
        A magnetic field. Must be a 3-component array or list.
        Can only be passed when using the Microwave channel ("mw_global").
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

    .. jupyter-execute::

        import qse
        qbits = qse.lattices.square(4.0, 3, 3)
        duration = 400
        amp = qse.Signal([1.01], duration)
        det = qse.Signal([0.12], duration)
        pcalc = qse.calc.Pulser(
            qbits=qbits,
            amplitude=amp,
            detuning=det)
        pcalc.build_sequence()
        pcalc.calculate()
        pcalc.get_spins()

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
        channel="rydberg_global",
        magnetic_field=None,
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
            qbits=qbits,
            label=label,
            is_calculator_available=CALCULATOR_AVAILABLE,
            installation_message=installation_message,
        )
        self.device = pulser.devices.MockDevice if device is None else device
        self.emulator = QutipEmulator if emulator is None else emulator
        self.wtimes = wtimes
        self._results = None
        self.channel = channel
        self.magnetic_field = magnetic_field

        if magnetic_field is not None and channel != "mw_global":
            raise Exception(
                "A magnetic field can only be passed when using the Microwave channel."
            )

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
    def register(self):
        return self._qbits.to_pulser()

    @property
    def sequence(self):
        return self._sequence

    @property
    def results(self) -> SimResult:
        """
        By simply overriding the property signature and calling super(), 
        Pylance instantly knows this branch returns a SimResult.
        """
        return super().results  # type: ignore

    def build_sequence(self):
        """
        Build the sequence of operations involving the qubit coordinates,
        amplitude pulse and detuning pulse.
        """
        self._sequence = pulser.Sequence(self.register, self.device)
        self._sequence.declare_channel("ch0", self.channel)
        if self.magnetic_field is not None:
            self._sequence.set_magnetic_field(*self.magnetic_field)
        self._sequence.add(
            pulser.Pulse(amplitude=self.amplitude, detuning=self.detuning, phase=0),
            "ch0",
        )
        if self.channel == "mw_global":
            basis = "XY"
        else:
            basis = "ground-rydberg"
        self._sequence.measure(basis)
        self._sim = self.emulator.from_sequence(self._sequence)

    @property
    def sim(self):
        return self._sim

    def calculate(self, progress=True) -> SimResult:
        """
        Run the calculation.
        """
        # we need to have/add an attribute to calc for device
        # check whether all the emulator `classes` have .from_sequence method.
        if self.wtimes:
            t1 = time()

        # Execute the pulser backend
        raw_results = self.sim.run(progress_bar=progress)

        # Extract metadata from self.sim
        metadata = extract_safe_metadata(self.sim, "pulser")
        # metadata = {
        #    "backend": "pulser",
        #    "basis": getattr(self.sim, "basis_name", "unknown"),
        #    "total_duration_ns": getattr(self.sim, "total_duration_ns", -1),
        #    "n_time_steps": len(self.sim.evaluation_times),
        #    "has_noise": hasattr(self.sim, "noise_model")
        #    and self.sim.noise_model is not None,
        #    "noise_model": getattr(self.sim, "noise_model", "unknown"),
        # }

        # bind local result to the global extractor functions using "partial"
        # This creates new functions/callables that already have results loaded
        # as their first argument.

        bind_statevector = partial(extract_pulser_statevector, raw_results)
        bind_counts = partial(extract_pulser_counts, raw_results)
        bind_expectation = partial(extract_pulser_expectation, raw_results)
        bind_states = partial(extract_pulser_states, raw_results)

        # In the qutip backend pulser uses the convention of 0 (1) being
        # the excited (ground) state. Hence we must reverse the state vector.
        # at the moment there does not seem an effective or better alternative than
        # self.spins = self.get_spins()

        if self.wtimes:
            t2 = time()
            exec_time = t2 - t1
            metadata["exec_time"] = exec_time
            print(f"time in compute and simulation = {exec_time} s.")

        # Return the clean SimResult
        # SimResult can now call bound_statevector() with zero arguments!
        self._results = SimResult(
            statevector_func=bind_statevector,
            counts_func=bind_counts,
            expectation_func=bind_expectation,
            states_generator=bind_states,
            metadata=metadata)
        return self.results



def _format_pulse(pulse):
    if pulse is None or isinstance(pulse, pulser.waveforms.Waveform):
        return pulse
    if isinstance(pulse, (Signal, Signals)):
        return pulse.to_pulser()

    raise Exception(
        "Pulses must be either `qse.Signal`, `qse.Signals` "
        "or `pulser.waveforms.Waveform`."
    )
