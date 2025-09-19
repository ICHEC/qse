import numpy as np
import pytest

import qse

pauli_x = np.array([[0.0, 1.0], [1.0, 0.0]])
pauli_y = np.array([[0.0, -1.0j], [1.0j, 0.0]])
pauli_z = np.array([[1.0, 0.0], [0.0, -1.0]])


def random_state(hdim):
    statevector = np.random.rand(hdim) + 1j * np.random.rand(hdim)
    statevector /= np.linalg.norm(statevector)
    return statevector


def kron_list(ops):
    x = ops[0]
    if len(ops) == 1:
        return x
    for op in ops[1:]:
        x = np.kron(x, op)
    return x


def expval(op, state):
    return (np.conj(state) @ (op @ state[:, None])).item()


@pytest.mark.parametrize("n_qubits, hsize", [(4, 2**4), (3, 2**3-1)])
def test_get_basis_shape(n_qubits, hsize):
    """Check that get_basis outputs a n_qubits * h_size shape"""
    ibasis = qse.magnetic.get_basis(n_qubits, hsize)
    assert ibasis.shape == (hsize, n_qubits)

@pytest.mark.parametrize("hdim", [1,2,3,4,5])
def test_spin_values_are_less_than_one(hdim):
    "Test that the absolute values of the computed spins is less than 1"
    statevector = random_state(2**hdim)
    basis = qse.magnetic.get_basis(hdim, 2**hdim)
    spins = qse.magnetic.get_spins(statevector, hdim, basis)
    assert np.all(np.abs(spins) <= 1)


def test_spin_two_qubits():
    """Check get_spins works for two qubits."""
    n = 2
    statevector = random_state(2**n)
    spin_qse = qse.magnetic.get_spins(statevector, n)

    op_q0 = [np.kron(p, np.eye(2)) for p in [pauli_x, pauli_y, pauli_z]]
    expval_q0 = [expval(op, statevector) for op in op_q0]
    op_q1 = [np.kron(np.eye(2), p) for p in [pauli_x, pauli_y, pauli_z]]
    expval_q1 = [expval(op, statevector) for op in op_q1]

    spin_np = np.stack([expval_q0, expval_q1])

    assert np.allclose(spin_qse, spin_np)


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
    psi = random_state(2**n)

    np_n_op = np.zeros(n, dtype=complex)
    for i in range(n):
        op_list = [np.eye(2)] * n
        op_list[i] = np.array([[0.0, 0.0], [0.0, 1.0]])
        np_n_op[i] = expval(kron_list(op_list), psi)

    n_op = qse.magnetic.get_number_operator(psi, n)
    assert np.allclose(n_op, np_n_op)
