import numpy as np

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
