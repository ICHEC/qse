import numpy as np
import pytest

import qse


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
    assert all(n_op[i] == float(basis_state[i]) for i in range(n))


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
    assert all(n_op[i] == np_n_op[i] for i in range(n))
