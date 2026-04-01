import numpy as np
import qutip as qp

from qse import Operator, Operators
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
        self.amplitude = amplitude
        self.detuning = detuning

        if self.amplitude.duration != self.detuning.duration:
            raise Exception("Amplitude and Detuning signals must have same duration.")

        self._build_hamiltonians()

    def calculate(self, e_ops=None):
        state = qp.tensor([qp.fock(2)] * self.qbits.nqbits)

        if e_ops is not None:
            exps = []

        for i in range(self.amplitude.duration):
            print(i)
            hamiltonian = self.get_hamiltonian(i)
            evo_step = qp.sesolve(hamiltonian, state, [0, 1e-3])
            state = evo_step.final_state
            if e_ops is not None:
                exps.append([qp.expect(op, state) for op in e_ops])

        result = {"state": state}
        if e_ops is not None:
            result["expectations" : np.array(exps)]

        return result

    def get_hamiltonian(self, i):
        return (
            0.5 * self.ham_amplitude * self.amplitude.values[i]
            - self.ham_detuning * self.detuning.values[i]
            + self.ham_int
        )

    def _build_hamiltonians(self):
        interaction = self.qbits.compute_interaction_hamiltonian(
            lambda d: self.c6 / d**6, "N"
        )
        self.ham_int = interaction.to_qutip()

        n = self.qbits.nqbits
        self.ham_amplitude = Operators(
            [Operator("X", i, n, 0.5) for i in range(n)]
        ).to_qutip()
        self.ham_detuning = Operators(
            [Operator("N", i, n) for i in range(n)]
        ).to_qutip()
