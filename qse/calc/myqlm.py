"""

Myqlm Calculator
----------------

This module Provides QSE calculator interface to MyQLM backend.

Ref: https://myqlm.github.io/

`myqlm` is a software middleware to potentially several quantum
hardware providers, such as Pasqal devices for analog quantum
computing.

It lets one offload an analog quantum workflow to myqlm middleware
which could be simulated ot real quantum device.

It has been tested in HPC setting for simulated backend, and
will be tested when we get access to a real Pasqal device.

"""

from functools import partial
from time import time
from typing import Any

import numpy as np

from qse.calc.calculator import Calculator
from qse.calc.results import SimResult

qat_available = False
qlmaas_available = False

try:
    import qat.qpus

    if hasattr(qat.qpus, "AnalogQPU"):
        AQPU_local = qat.qpus.AnalogQPU
        qat_available = True
    else:
        AQPU_local = None
except ModuleNotFoundError:
    qat_available = False

try:
    import qlmaas.qpus

    if hasattr(qlmaas.qpus, "AnalogQPU"):
        AQPU_remote = qlmaas.qpus.AnalogQPU
        qlmaas_available = True
    else:
        AQPU_remote = None
except (ModuleNotFoundError, ImportError):
    qlmaas_available = False

if qat_available:
    AQPU = AQPU_local
else:
    if qlmaas_available:
        AQPU = AQPU_remote
    else:
        AQPU = None

try:
    import qat

    CALCULATOR_AVAILABLE = True
except ImportError:
    CALCULATOR_AVAILABLE = False


default_params = {
    "rydberg": {
        "C6": 5420158.53,  # unit (rad/µs)(µm)**6
        "min_omega": None,
        "max_omega": 12.57,  # rad/µs from pulser
        "min_delta": None,
        "max_delta": 125.7,  # rad/µs from pulser
        "min_atom_distance": 4,  # µm from pulser
        "max_duration": 4000,  # ns from pulser (max sequence duration),
    },
    "ssh": {},
    "fermion": {},
}


def extract_myqlm_statevector(raw_results: Any, num_qubits: int) -> np.ndarray:
    """
    Reconstructs the full 2^N dense statevector from myQLM's truncated iterable.
    Guarantees the output array size matches the Hilbert space dimension.
    """
    # Guard clause: No amplitudes returned (e.g., in purely noisy runs)
    if not hasattr(raw_results[0], "amplitude"):
        return None 
    hsize = 2**num_qubits
    
    # 1. Extract the surviving basis integers and their complex amplitudes
    basis = np.fromiter((s.state.int for s in raw_results), dtype=int)
    amplitudes = np.fromiter((s.amplitude for s in raw_results), dtype=complex)
    
    # 2. Reconstruct the full Hilbert space (your zero-padding logic)
    if amplitudes.shape[0] < hsize:
        statevector = np.zeros(hsize, dtype=complex)
        # Advanced indexing drops the surviving amplitudes into their exact integer positions
        statevector[basis] = amplitudes
        return statevector
        
    return amplitudes


def extract_myqlm_counts(raw_results: Any, shots: int) -> dict:
    """
    Map the bitstring basis states to their probabilities.
    """
    counts = {}
    for sample in raw_results:
        # Convert the integer state to a zero-padded bitstring
        # bitstring = format(sample.state.int, f'0{num_qubits}b')
        bit_int = sample.state.int
        counts[bit_int] = int(sample.probability * shots)
    return counts


def extract_myqlm_expectation(raw_results) -> float:
    """
    Extracts the expectation value if the job was an observable measurement.
    """
    # If the job was an observable (not a standard statevector simulation),
    # myQLM stores the expectation result in the .value attribute.
    if hasattr(raw_results, 'value') and raw_results.value is not None:
        return raw_results.value
        
    # If it was a standard sampling/statevector job, expectation is N/A here.
    return None


def extract_myqlm_states(raw_results: Any) -> Any:
    """
    Returns the basis states as a list of bitstrings.
    To match the zero-padded statevector, we generate the full 2^N Hilbert space basis.
    """
    if hasattr(raw_results, 'states'):
        states = raw_results.states
    else:
        states = None
    # Generating the full basis ensures it perfectly aligns with the zero-padded statevector
    # [format(i, f'0{num_qubits}b') for i in range(2**num_qubits)]
    return states


class Myqlm(Calculator):
    """
    QSE-Calculator for MyQLM.

    Parameters
    ----------
    qbits: qse.Qbits
        The qbits object.
    amplitude : qse.Signal
        The amplitude pulse.
    detuning : qse.Signal
        The detuning pulse.
    qpu : qat.qpus
        The Quantum Processing Unit for executing the job.
    label: str
        The label.
        Defaults to "pulser-run".
    wtimes: bool
        Whether to print the times.
        Defaults to True.
    """

    def __init__(
        self,
        qbits=None,
        amplitude=None,
        detuning=None,
        qpu=None,
        label="myqlm-run",
        wtimes=True,
        shots=1024
    ):
        installation_message = (
            "myQLM is not installed. To install, "
            "see https://myqlm.github.io/01_getting_started/:myqlm:01_install.html."
        )

        super().__init__(
            qbits=qbits,
            label=label,
            is_calculator_available=CALCULATOR_AVAILABLE,
            installation_message=installation_message,
        )

        self.qpu = AQPU() if qpu is None else qpu
        self.label = label
        self.wtimes = wtimes
        self.shots = shots
        self._results = None
        self.system = "rydberg"

        self.params = dict(default_params[self.system])

        self.amplitude = amplitude
        self.detuning = detuning

        self.C6 = self.params["C6"]

    @property
    def Hamiltonian(self):
        return self._get_hamiltonian()

    def _get_hamiltonian(self):
        if self.system == "rydberg":
            ham = self._generate_rydberg_hamiltonian()
        else:
            ham = None
        return ham

    def _generate_rydberg_hamiltonian(self):
        nqbits = self.qbits.nqbits
        # > add checks for ensuring whether amplitudes/detunings etc
        # > are within allowed limits compatible with pulser virtual device
        H_amplitude = qat.core.Observable(
            nqbits, pauli_terms=[qat.core.Term(0.5, "X", [i]) for i in range(nqbits)]
        )
        amplitude = _waveform(self.amplitude)

        H_detuning = qat.core.Observable(
            nqbits, pauli_terms=[qat.core.Term(0.5, "Z", [i]) for i in range(nqbits)]
        )
        detuning = _waveform(self.detuning)

        rij = self.qbits.get_all_distances()
        H_interact = 0
        for i in range(nqbits):
            for j in range(i + 1, nqbits):
                H_interact += (
                    (self.C6 / rij[i, j] ** 6) * _occ_op(nqbits, i) * _occ_op(nqbits, j)
                )

        return [
            (amplitude, H_amplitude),
            (detuning, H_detuning),
            (1, H_interact),
        ]

    @property
    def results(self) -> SimResult:
        """
        By simply overriding the property signature and calling super(),
        Pylance instantly knows this branch returns a SimResult.
        """
        return super().results  # type: ignore
    
    def calculate(self):
        """
        Run the calculation.
        """
        if self.wtimes:
            t1 = time()

        tmax = (
            max(self.amplitude.duration, self.detuning.duration) / 1000
        )  # convert to microseconds.
        self.schedule = qat.core.Schedule(drive=self.Hamiltonian, tmax=tmax)
        self.job = self.schedule.to_job()
        self.async_result = self.qpu.submit(self.job)
        raw_results = self.async_result.join()

        # get the number of qubits
        num_qubits = raw_results.nbqbits or len(self.qbits)

        # Bind the parsers

        bind_statevector = partial(extract_myqlm_statevector, raw_results, num_qubits)
        bind_counts = partial(extract_myqlm_counts, raw_results)
        bind_expectation = partial(extract_myqlm_expectation, raw_results)
        bind_states = partial(extract_myqlm_states, raw_results)        

        # 2. Build the metadata from the wrapper attributes
        metadata = {
            "backend": "myqlm",
            "num_qubits": num_qubits,
            "has_statevector": getattr(raw_results, 'has_statevector', False),
            "myqlm_metadata": getattr(raw_results, 'meta_data', {}),
            "raw_results": raw_results
        }

        # 3. Instantiate the Result
        self._results = SimResult(
            statevector_func=bind_statevector,
            counts_func=bind_counts,
            expectation_func=bind_expectation,
            states_generator=bind_states,
            metadata=metadata,
        )

        # self.spins = self.get_spins()
        ## self.sij = self.get_sij()
        if self.wtimes:
            t2 = time()
            print(f"time in compute and simulation = {t2 - t1} s.")

        return self.results


def _waveform(pulse):
    vi = pulse.values
    tmax = pulse.duration / 1000  # convert to microseconds.
    ti = np.linspace(0, tmax, vi.shape[0])
    vi_m = np.diff(vi)
    ti_m = np.diff(ti)
    vi_p = vi[1:] + vi[:-1]
    ti_p = ti[1:] + ti[:-1]
    a = vi_m / ti_m
    b = 0.5 * (vi_p - a * ti_p)
    arith_expr = 0
    t_var = qat.core.Variable("t")
    for i, (ai, bi) in enumerate(zip(a, b)):
        # Create ax + b by calculating the slope and the offset
        respective_line = ai * t_var + bi
        arith_expr += (
            qat.core.variables.heaviside(t_var, ti[i], ti[i + 1]) * respective_line
        )
    return arith_expr


def _occ_op(nqbits, qi):
    ti = qat.core.Term(1.0, "Z", [qi])
    return (1 - qat.core.Observable(nqbits, pauli_terms=[ti])) / 2
