"""
This module Provides QSE calculator interface to MyQLM. It is derived from
the ASE's calculator, though at the end it may look very different from
ASE calculator.
https://myqlm.github.io/
"""

from time import time

import numpy as np

from qse.calc.calculator import Calculator

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
    "rydberg_global": {
        "C6": 5420158.53,  # unit (rad/µs)(µm)**6
        "min_omega": None,
        "max_omega": 12.57,  # rad/µs from pulser
        "min_delta": None,
        "max_delta": 125.7,  # rad/µs from pulser
        "min_atom_distance": 4,  # µm from pulser
        "max_duration": 4000,  # ns from pulser (max sequence duration),
    },
    "mw_global": {
        "C3": 3700,  # unit (rad/µs)(µm)**3
    },
}


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
    channel : str
        Which channel to use. For example "rydberg_global" for Rydberg or
        "mw_global" for microwave.
        Defaults to "rydberg_global".
    magnetic_field : np.ndarray | list
        A magnetic field. Must be a 3-component array or list.
        Can only be passed when using the Microwave channel ("mw_global").
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
        channel="rydberg_global",
        magnetic_field=None,
        qpu=None,
        label="myqlm-run",
        wtimes=True,
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
        self.wtimes = wtimes
        self.results = None

        self.channel = channel
        self.magnetic_field = magnetic_field
        self.params = dict(default_params[self.channel])

        self.amplitude = amplitude
        self.detuning = detuning

    def get_hamiltonian(self):
        if self.system == "rydberg_global":
            return self._generate_rydberg_hamiltonian(
                self.qbits, self.amplitude, self.detuning, self.params["C6"]
            )
        elif self.system == "mw_global":
            return self._generate_microwave_hamiltonian(
                self.qbits,
                self.amplitude,
                self.detuning,
                self.params["C3"],
                self.magnetic_field,
            )
        return None

    def calculate(self):
        """
        Run the calculation.
        """
        if self.wtimes:
            t1 = time()

        tmax = (
            max(self.amplitude.duration, self.detuning.duration) / 1000
        )  # convert to microseconds.
        self.schedule = qat.core.Schedule(drive=self.get_hamiltonian(), tmax=tmax)
        self.job = self.schedule.to_job()
        self.async_result = self.qpu.submit(self.job)
        self.results = self.async_result.join()

        self.probabities = np.fromiter(
            (s.probability for s in self.results), dtype=float
        )
        self.basis = np.fromiter((s.state.int for s in self.results), dtype=int)

        # don't know why result.statevector is Nonetype, fill it with state
        if hasattr(self.results[0], "amplitude"):
            statevector = np.fromiter(
                (s.amplitude for s in self.results), dtype=complex
            )
            N = len(self.qbits)
            hsize = 2**N
            if statevector.shape[0] < hsize:
                coeff0 = np.zeros(hsize, dtype=complex)
                coeff0[self.basis] = statevector
                statevector = coeff0
            self.statevector = statevector
        self.spins = self.get_spins()
        # self.sij = self.get_sij()
        if self.wtimes:
            t2 = time()
            print(f"time in compute and simulation = {t2 - t1} s.")


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


def _plus_op(nqbits, qi):
    return (
        qat.core.Observable(nqbits, pauli_terms=[qat.core.Term(1.0, "X", [qi])])
        - 1j * qat.core.Observable(nqbits, pauli_terms=[qat.core.Term(1.0, "Y", [qi])])
    ) / 2


def _minus_op(nqbits, qi):
    return (
        qat.core.Observable(nqbits, pauli_terms=[qat.core.Term(1.0, "X", [qi])])
        + 1j * qat.core.Observable(nqbits, pauli_terms=[qat.core.Term(1.0, "Y", [qi])])
    ) / 2


def _hopping_op(nqbits, qi, qj):
    return _plus_op(nqbits, qi) * _minus_op(nqbits, qj) + _plus_op(
        nqbits, qj
    ) * _minus_op(nqbits, qi)


def _generate_rydberg_hamiltonian(qbits, amplitude, detuning, c6):
    nqbits = qbits.nqbits
    H_amplitude = qat.core.Observable(
        nqbits, pauli_terms=[qat.core.Term(0.5, "X", [i]) for i in range(nqbits)]
    )
    amplitude = _waveform(amplitude)

    H_detuning = qat.core.Observable(
        nqbits, pauli_terms=[qat.core.Term(0.5, "Z", [i]) for i in range(nqbits)]
    )
    detuning = _waveform(detuning)

    rij = qbits.get_all_distances()
    H_interact = 0
    for i in range(nqbits):
        for j in range(i + 1, nqbits):
            H_interact += (
                (c6 / rij[i, j] ** 6) * _occ_op(nqbits, i) * _occ_op(nqbits, j)
            )

    return [
        (amplitude, H_amplitude),
        (detuning, H_detuning),
        (1, H_interact),
    ]


def _generate_microwave_hamiltonian(qbits, amplitude, detuning, c3, magnetic_field):
    nqbits = qbits.nqbits
    H_amplitude = qat.core.Observable(
        nqbits, pauli_terms=[qat.core.Term(0.5, "X", [i]) for i in range(nqbits)]
    )
    amplitude = _waveform(amplitude)

    H_detuning = qat.core.Observable(
        nqbits, pauli_terms=[qat.core.Term(0.5, "Z", [i]) for i in range(nqbits)]
    )
    detuning = _waveform(detuning)

    rij = qbits.get_all_distances()

    if magnetic_field is None:
        H_interact = 0
        for i in range(nqbits):
            for j in range(i + 1, nqbits):
                H_interact += (c3 / rij[i, j] ** 3) * _hopping_op(nqbits, i, j)
    else:
        H_interact = 0
        magnetic_field = np.array(magnetic_field)
        magnetic_field /= np.linalg.norm(magnetic_field)  # normalize
        for i in range(nqbits):
            for j in range(i + 1, nqbits):
                qubit_vec = qbits.positions[i] - qbits.positions[j]
                qubit_vec /= np.linalg.norm(qubit_vec)
                m_term = 1 - 3 * np.dot(magnetic_field, qubit_vec) ** 2
                H_interact += (m_term * c3 / rij[i, j] ** 3) * _hopping_op(nqbits, i, j)

    return [
        (amplitude, H_amplitude),
        (detuning, H_detuning),
        (1, H_interact),
    ]
