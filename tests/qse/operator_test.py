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
    "indices, operator, expected_str",
    [
        (0, "X", "XI"),
        (0, ["X"], "XI"),
        (1, "Y", "IY"),
        ([0, 1], "Z", "ZZ"),
        ([0, 1], ["X", "Y"], "XY"),
        ([1, 0], ["X", "Y"], "YX"),
    ],
)
def test_operator_2qubits(indices, operator, expected_str):
    """Test creating an Operator with 2 qubits."""
    op = qse.Operator(operator, indices, nqbits=2)

    assert isinstance(op, qse.Operator)
    assert op.to_str() == expected_str
    assert np.isclose(op.coef, 1.0)

    assert np.allclose(
        np.kron(pauli_dict[expected_str[0]], pauli_dict[expected_str[1]]),
        op.to_qutip().full(),
    )


@pytest.mark.parametrize(
    "indices, operator, expected_str",
    [
        (0, "X", "XIII"),
        (1, ["X"], "IXII"),
        ([0, 1], "Z", "ZZII"),
        ([1, 2], ["Z", "Y"], "IZYI"),
    ],
)
def test_operator_4qubits(indices, operator, expected_str):
    """Test creating an Operator with 4 qubits."""
    op = qse.Operator(operator, indices, nqbits=4)

    assert isinstance(op, qse.Operator)
    assert op.to_str() == expected_str
    assert np.isclose(op.coef, 1.0)


@pytest.mark.parametrize("op1", ["X", "Y"])
@pytest.mark.parametrize("op2", ["Z", "Y"])
@pytest.mark.parametrize("coef", [0.1, -0.45, 1j + 53])
def test_operator_4qubits_coef(op1, op2, coef):
    """Test creating an Operator with 4 qubits and a coefficient."""
    op = qse.Operator(operator=[op1, op2], indices=[1, 2], nqbits=4, coef=coef)

    assert isinstance(op, qse.Operator)
    assert op.to_str() == "I" + op1 + op2 + "I"
    assert np.isclose(op.coef, coef)

    op_np = coef * np.kron(
        np.kron(np.kron(pauli_dict["I"], pauli_dict[op1]), pauli_dict[op2]),
        pauli_dict["I"],
    )
    assert np.allclose(op_np, op.to_qutip().full())


@pytest.mark.parametrize(
    "indices, operator",
    [
        (0, "x"),
        (0, ["X", "Y"]),
        (0, "XY"),
    ],
)
def test_operator_fail(indices, operator):
    """Test creating an Operator with 4 qubits."""
    with pytest.raises(Exception):
        qse.Operator(operator, indices, nqubits=4)


@pytest.mark.parametrize("qubits", [-1, 4, [-1, 2], [3, 4]])
def test_operator_fail_qubit_index(qubits):
    """Test the class raises an error for a qubit outside the index range."""
    with pytest.raises(Exception):
        qse.Operator("Z", qubits, nqubits=4)


def test_operator_mul():
    """Test multiplying an Operator by a scalar."""
    op = qse.Operator("X", 0, nqbits=2, coef=1.0)

    op_scaled = op * 3.5

    assert isinstance(op_scaled, qse.Operator)
    assert op_scaled.to_str() == "XI"
    assert np.isclose(op_scaled.coef, 3.5)

    op *= 2.0
    assert np.isclose(op.coef, 2.0)


def test_operators_init():
    """Test initializing empty Operators."""
    ops = qse.Operators()

    assert isinstance(ops, qse.Operators)
    assert ops.nterms == 0


def test_operators_creation_addition():
    """Test creating and adding Operators."""
    nqbits = 4
    ops = qse.Operators([qse.Operator("X", i, nqbits) for i in range(nqbits)])
    assert isinstance(ops, qse.Operators)
    assert ops.nterms == nqbits
    assert ops.nqbits == nqbits

    ops += qse.Operators([qse.Operator("Y", i, nqbits) for i in range(nqbits)])
    assert isinstance(ops, qse.Operators)
    assert ops.nterms == nqbits * 2
    assert ops.nqbits == nqbits

    for i in range(2):
        ops += qse.Operator("Z", i, nqbits)
        assert isinstance(ops, qse.Operators)
        assert ops.nterms == nqbits * 2 + i + 1
        assert ops.nqbits == nqbits


def test_operators_qutip():
    """Test to_qutip method."""
    nqbits = 3
    op1 = qse.Operator("X", 0, nqbits, 9.242)
    op2 = qse.Operator(["X", "Y"], [1, 2], nqbits, 1.0)
    ops = qse.Operators([op1, op2])

    op_sum = op1.to_qutip() + op2.to_qutip()
    assert np.allclose(op_sum.full(), ops.to_qutip().full())


def test_operators_mul():
    """Test multiplying Operators by a scalar."""
    nqbits = 4
    ops = qse.Operators(
        [qse.Operator("X", i, nqbits=nqbits + 5) for i in range(nqbits)]
    )
    assert all(np.isclose(op.coef, 1.0) for op in ops)

    ops_scaled = ops * 2.5
    assert all(np.isclose(op.coef, 2.5) for op in ops_scaled)

    ops *= 3.5
    assert all(np.isclose(op.coef, 3.5) for op in ops)
