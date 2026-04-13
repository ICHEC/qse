import numpy as np
import qutip as qp

from qse import Operator, Operators, Signals
from qse.calc.calculator import Calculator

c6 = 5420158.53  # (rad/µs)(µm)**6


class Qutip(Calculator):
    c6 = c6  # (rad/µs)(µm)**6

    def __init__(self, qbits, amplitude, detuning):
        super().__init__(
            qbits=qbits,
            label="qutip",
            is_calculator_available=True,
            installation_message="",
        )

        if not isinstance(amplitude, Signals):
            amplitude = Signals([amplitude])

        if not isinstance(detuning, Signals):
            detuning = Signals([detuning])

        _check_pulses(amplitude, detuning)

        self.amplitude = amplitude
        self.detuning = detuning

        self._build_hamiltonian_operators()

    def calculate(self, e_ops=None):
        state = qp.tensor([qp.fock(2)] * self.qbits.nqbits)

        if e_ops is not None:
            exps = []

        for amp_s, det_s in zip(self.amplitude, self.detuning):
            delta_t = _nano_to_micro(amp_s.time_per_value())
            for amp, det in zip(amp_s.values, det_s.values):
                hamiltonian = self.get_hamiltonian(amp, det)
                evo_step = qp.sesolve(hamiltonian, state, [0, delta_t])
                state = evo_step.final_state
                if e_ops is not None:
                    exps.append([qp.expect(op, state) for op in e_ops])

        result = {"state": state.full()}
        if e_ops is not None:
            result["expectations" : np.array(exps)]

        return result

    def get_hamiltonian(self, amplitude, detuning):
        return (
            self.ham_amplitude * amplitude
            + self.ham_detuning * detuning
            + self.ham_int
        )

    def _build_hamiltonian_operators(self):
        interaction = self.qbits.compute_interaction_hamiltonian(
            lambda d: self.c6 / d**6, "N"
        )
        self.ham_int = interaction.to_qutip()

        n = self.qbits.nqbits
        self.ham_amplitude = Operators(
            [Operator("X", i, n, 0.5) for i in range(n)]
        ).to_qutip()
        self.ham_detuning = Operators(
            [Operator("N", i, n, 0.5) for i in range(n)]
        ).to_qutip()


def _check_pulses(amplitude, detuning):
    if amplitude.duration != detuning.duration:
        raise Exception("The amplitude and detuning must have the same duration.")

    if len(amplitude) != len(detuning):
        raise Exception(
            "The amplitude and detuning must contain the same number of signals."
        )

    for a, d in zip(amplitude, detuning):
        if len(a) != len(d):
            raise Exception(
                "The amplitude and detuning must have the same amount of values."
            )


def _nano_to_micro(t):
    return t / 1000
