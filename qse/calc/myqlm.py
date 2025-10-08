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
    "rydberg": {
        "amplitude": None,
        "detuning": None,
        "C6": 5420158.53,  # unit (rad/µs)(µm)**6
        "min_omega": None,
        "max_omega": 12.57,  # rad/µs from pulser
        "min_delta": None,
        "max_delta": 125.7,  # rad/µs from pulser
        "min_atom_distance": 4,  # µm from pulser
        "max_duration": 4000,  # ns from pulser (max sequence duration),
        "default_points": 6,
    },
    "ssh": {},
    "fermion": {},
}


# for rydberg system we need additional parameters: mix,max of
# amplitude and detuning to set based on device.
class Myqlm(Calculator):
    """
    QSE-Calculator for MyQLM

    Parameters
    ----------
    qbits
        ...
    amplitude
        ...
    detuning
        ...
    qpu
        ...
    system
        ...
    label
        ...
    wtimes
        ...
    """

    implemented_properties = ["energy", "state", "fidality"]
    default_parameters = dict(label="q", qbits=None)

    def __init__(
        self,
        qbits=None,
        amplitude=None,
        detuning=None,
        qpu=None,
        system="rydberg",
        label="myqlm-run",
        wtimes=True,
    ):
        installation_message = (
            "myQLM is not installed. To install, "
            "see https://myqlm.github.io/01_getting_started/:myqlm:01_install.html."
        )

        super().__init__(
            CALCULATOR_AVAILABLE, installation_message, label=label, qbits=qbits
        )

        self.qpu = AQPU if qpu is None else qpu
        self.label = label
        self.wtimes = wtimes
        self.results = None
        self.system = system
        self.params = dict(default_params[self.system])

        self.duration = amplitude.duration / 1000  # convert to µs

        self.amp = (
            amplitude
            if amplitude is not None
            else np.zeros(self.params["default_points"])
        )
        self.amplitude = self._waveform(self.amp.values, tmax=self.duration)
        self.det = (
            detuning
            if detuning is not None
            else np.zeros(self.params["default_points"])
        )
        self.detuning = self._waveform(self.det.values, tmax=self.duration)

        self.C6 = self.params["C6"]
        self.qpu = None
        self.spins = None
        self.sij = None

    def _occ_op(self, nqbits, qi):
        ti = qat.core.Term(1.0, "Z", [qi])
        return (1 - qat.core.Observable(nqbits, pauli_terms=[ti])) / 2

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
        rij = self.qbits.get_all_distances()
        nqbits = self.qbits.nqbits
        # > add checks for ensuring whether amplitudes/detunings etc
        # > are within allowed limits compatible with pulser virtual device
        amplitude = self.amplitude
        detuning = self.detuning

        H1_terms = [qat.core.Term(0.5, "X", [i]) for i in range(nqbits)]
        H2_terms = [qat.core.Term(0.5, "Z", [i]) for i in range(nqbits)]
        H_amplitude = qat.core.Observable(nqbits, pauli_terms=H1_terms)
        H_detuning = qat.core.Observable(nqbits, pauli_terms=H2_terms)
        H_interact = 0
        for i in range(nqbits):
            for j in range(i + 1, nqbits):
                H_interact += (
                    (self.C6 / rij[i, j] ** 6)
                    * self._occ_op(nqbits, i)
                    * self._occ_op(nqbits, j)
                )

        return [
            (amplitude, H_amplitude),
            (detuning, H_detuning),
            (1, H_interact),
        ]

    def _waveform(self, vi, tmax):
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

    def calculate(self, qbits=None, properties=..., system_changes=...):
        """
        _summary_

        Parameters
        ----------
        qbits : _type_, optional
            _description_, by default None
        properties : _type_, optional
            _description_, by default ...
        system_changes : _type_, optional
            _description_, by default ...
        """
        if self.wtimes:
            t1 = time()

        self.schedule = qat.core.Schedule(drive=self.Hamiltonian, tmax=self.duration)
        self.job = self.schedule.to_job()
        self.qpu = AQPU()
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
            )[
                ::-1
            ]  # myqlm uses opposite statevector convention to normal.
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
