import numpy as np
import pytest

import qse


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
