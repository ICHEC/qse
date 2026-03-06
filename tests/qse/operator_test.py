import numpy as np
import pytest

import qse

pauli_dict = {
    "I": np.eye(2),
    "X": np.array([[0.0, 1.0], [1.0, 0.0]]),
    "Y": np.array([[0.0, -1.0j], [1.0j, 0.0]]),
    "Z": np.array([[1.0, 0.0], [0.0, -1.0]]),
}


@pytest.mark.parametrize(
    "qubits, operator, expected_str",
    [
        (0, "X", "XI"),
        (0, ["X"], "XI"),
        (1, "Y", "IY"),
        ([0, 1], "Z", "ZZ"),
        ([0, 1], ["X", "Y"], "XY"),
        ([1, 0], ["X", "Y"], "YX"),
    ],
)
def test_operator_2qubits(qubits, operator, expected_str):
    """Test creating an Operator with 2 qubits."""
    op = qse.Operator(operator, qubits, nqubits=2)

    assert isinstance(op, qse.Operator)
    assert op.to_str() == expected_str
    assert np.isclose(op.coef, 1.0)

    assert np.allclose(
        np.kron(pauli_dict[expected_str[0]], pauli_dict[expected_str[1]]),
        op.to_qutip().full(),
    )


@pytest.mark.parametrize(
    "qubits, operator, expected_str",
    [
        (0, "X", "XIII"),
        (1, ["X"], "IXII"),
        ([0, 1], "Z", "ZZII"),
        ([1, 2], ["Z", "Y"], "IZYI"),
    ],
)
def test_operator_4qubits(qubits, operator, expected_str):
    """Test creating an Operator with 4 qubits."""
    op = qse.Operator(operator, qubits, nqubits=4)

    assert isinstance(op, qse.Operator)
    assert op.to_str() == expected_str
    assert np.isclose(op.coef, 1.0)


@pytest.mark.parametrize("op1", ["X", "Y"])
@pytest.mark.parametrize("op2", ["Z", "Y"])
@pytest.mark.parametrize("coef", [0.1, -0.45, 1j + 53])
def test_operator_4qubits_coef(op1, op2, coef):
    """Test creating an Operator with 4 qubits and a coefficient."""
    op = qse.Operator(operator=[op1, op2], qubits=[1, 2], nqubits=4, coef=coef)

    assert isinstance(op, qse.Operator)
    assert op.to_str() == "I" + op1 + op2 + "I"
    assert np.isclose(op.coef, coef)

    op_np = coef * np.kron(
        np.kron(np.kron(pauli_dict["I"], pauli_dict[op1]), pauli_dict[op2]),
        pauli_dict["I"],
    )
    assert np.allclose(op_np, op.to_qutip().full())


@pytest.mark.parametrize(
    "qubits, operator",
    [
        (0, "x"),
        (0, ["X", "Y"]),
        (0, "XY"),
    ],
)
def test_operator_fail(qubits, operator):
    """Test creating an Operator with 4 qubits."""
    with pytest.raises(Exception):
        qse.Operator(operator, qubits, nqubits=4)
