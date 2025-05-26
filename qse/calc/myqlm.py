"""
This module Provides QSE calculator interface to MyQLM. It is derived from
the ASE's calculator, though at the end it may look very different from
ASE calculator.
https://myqlm.github.io/
"""

import numpy as np

from qse import Signal
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
except ModuleNotFoundError:
    qlmaas_available = False

if qat_available:
    AQPU = AQPU_local
else:
    if qlmaas_available:
        AQPU = AQPU_remote
    else:
        AQPU = None

if qat_available:
    print("qat is available")
else:
    print("qat is not available")
if qlmaas_available:
    print("qlmaas is available")
else:
    print("qlmaas is not available")
print(AQPU)

from time import time

# analogQPU imported based on what's available
import qat

import qse.magnetic as magnetic

# from qat.core.variables import Variable, heaviside


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
        "default_duration": 0.5,
        "default_points": 6,
    },
    "ssh": {},
    "fermion": {},
}


# for rydberg system we need additional parameters: mix,max of amplitude and detuning to set based on device.
class Myqlm(Calculator):
    """QSE-Calculator for MyQLM"""

    implemented_properties = ["energy", "state", "fidality"]
    default_parameters = dict(label="q", qbits=None)

    def __init__(
        self,
        qbits=None,
        amplitude=None,
        detuning=None,
        duration=None,
        qpu=None,
        analog=True,
        system="rydberg",
        label="myqlm-run",
        wtimes=True,
    ):
        super().__init__(label=label, qbits=qbits)
        self.qpu = AQPU if qpu is None else qpu
        self.label = label
        self.wtimes = wtimes
        self.results = None
        self.system = system
        self.params = dict(default_params[self.system])
        # if any of the default parameters are passed
        # as arguments, use them to update self.params
        self.duration = (
            duration if duration is not None else self.params["default_duration"]
        )
        self.amp = (
            amplitude
            if amplitude is not None
            else np.zeros(self.params["default_points"])
        )
        self.amplitude = self._waveform(self.amp, tmax=self.duration)
        self.det = (
            detuning
            if detuning is not None
            else np.zeros(self.params["default_points"])
        )
        self.detuning = self._waveform(self.det, tmax=self.duration)
        # self.duration = len(self.amplitude) # needs to change for time independent problems
        self.C6 = self.params["C6"]
        #
        self.qpu = None

    def _occ_op(self, nqbits, qi):
        ti = qat.core.Term(1.0, "Z", [qi])
        return (1 + qat.core.Observable(nqbits, pauli_terms=[ti])) / 2

    @property
    def Hamiltonian(self):
        return self._get_hamiltonian()

    def _get_hamiltonian(self):
        if self.system == "rydberg":
            ham = self._generate_rydberg_hamiltonian()
        else:
            ham = None
        return ham

    #

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
        Hamiltonian = [  # qat.core.Observable(
            (amplitude, H_amplitude),
            (detuning, H_detuning),
            (1, H_interact),
        ]
        return Hamiltonian

    #

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
        # return super().calculate(qbits, properties, system_changes)
        # self.Hamiltonian = self._get_hamiltonian()
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
            )
            N = len(self.qbits)
            hsize = 2**N
            ibasis = magnetic.get_basis(hsize=hsize, N=N)
            if statevector.shape[0] < hsize:
                coeff0 = np.zeros(hsize, dtype=complex)
                coeff0[self.basis] = statevector
                statevector = coeff0
            self.statevector = statevector
            del ibasis
        self.spins = self.get_spins()
        self.sij = self.get_sij()
        if self.wtimes:
            t2 = time()
            print(f"time in compute and simulation = {t2 - t1} s.")

    #

    def get_spins(self):
        """Get spin expectation values
        If the hamiltonian isn't simulated, it triggers simulation first.

        Returns:
            np.ndarray: Array of Nx3 containing spin expectation values.
        See :py.func: `qse.magnetic.get_spins` for more details.
        """
        if self.results is None:
            self.calculate()
        #
        nqbits = len(self.qbits)
        ibasis = magnetic.get_basis(2**nqbits, nqbits)
        si = magnetic.get_spins(self.statevector, ibasis, nqbits)
        return si

    def get_sij(self):
        r"""Get spin correlation s_ij
        If the hamiltonian isn't simulated, it triggers simulation first.

        Returns:
            np.ndarray: Array of NxN shape containing spin correlations.
        See :py.func: `qse.magnetic.get_sij` for more details.
        """
        if self.results is None:
            self.calculate()
        #
        nqbits = len(self.qbits)
        ibasis = magnetic.get_basis(2**nqbits, nqbits)
        sij = magnetic.get_sisj(self.statevector, ibasis, nqbits)
        return sij

    def structure_factor_from_sij(self, L1: int, L2: int, L3: int):
        r"""Get the structure factor

        Args:
            L1 (int): Extent of lattice in x direction
            L2 (int): Extent of lattice in y direction
            L3 (int): Extent of lattice in z direction

        Returns:
            np.ndarray: Array containing the structure factor
        See :py.func: `qse.magnetic.structure_factor_from_sij` for more details.
        """
        struc_fac = magnetic.structure_factor_from_sij(L1, L2, L3, self.qbits, self.sij)
        return struc_fac
