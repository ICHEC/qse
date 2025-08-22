import numpy as np
import pytest

import qse

pauli_x = np.array([[0.0, 1.0], [1.0, 0.0]])
pauli_y = np.array([[0.0, -1.0j], [1.0j, 0.0]])
pauli_z = np.array([[1.0, 0.0], [0.0, -1.0]])


def test_spin_two_qubits():
    """Check get_spins works for two qubits."""
    n = 2
    hdim = 2**n
    statevector = np.random.rand(hdim) + 1j * np.random.rand(hdim)
    statevector /= np.linalg.norm(statevector)

    def expval(ops, state):
        return [(np.conj(state) @ (op @ state[:, None])).item() for op in ops]

    spin_qse = qse.magnetic.get_spins(statevector, n)

    op_q0 = [np.kron(p, np.eye(2)) for p in [pauli_x, pauli_y, pauli_z]]
    op_q1 = [np.kron(np.eye(2), p) for p in [pauli_x, pauli_y, pauli_z]]
    spin_np = np.stack([expval(op_q0, statevector), expval(op_q1, statevector)])

    assert np.allclose(spin_qse, spin_np)


def kron_list(ops):
    x = ops[0]
    if len(ops) == 1:
        return x
    for op in ops[1:]:
        x = np.kron(x, op)
    return x


@pytest.mark.parametrize("basis_state", ["1", "00", "001", "101011"])
def test_number_operator_basis(basis_state):
    """Check the number operator for basis states."""
    n = len(basis_state)
    psi = np.zeros(2**n)
    psi[int(basis_state, 2)] = 1.0
    n_op = qse.magnetic.get_number_operator(psi, n)
    assert np.allclose(n_op, [float(basis_state[i]) for i in range(n)])


@pytest.mark.parametrize("n", [1, 2, 3])
def test_number_operator_random(n):
    """Check the number operator for random states."""
    psi = np.random.rand(2**n) + 1j * np.random.rand(2**n)
    psi /= np.linalg.norm(psi)

    np_n_op = np.zeros(n, dtype=complex)
    for i in range(n):
        op_list = [np.eye(2)] * n
        op_list[i] = np.array([[0.0, 0.0], [0.0, 1.0]])
        np_n_op[i] = (np.conj(psi) @ (kron_list(op_list) @ psi[:, None])).item()

    n_op = qse.magnetic.get_number_operator(psi, n)
    assert np.allclose(n_op, np_n_op)
