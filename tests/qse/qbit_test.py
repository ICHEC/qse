import numpy as np
import pytest

import qse


def test_default():
    """Check default behaviour."""
    qbit = qse.Qbit()
    assert isinstance(qbit, qse.Qbit)
    assert qbit.label == "X"
    assert np.allclose(qbit.position, np.zeros(3))


@pytest.mark.parametrize("label", ["a" * 3, "123a", None, 14])
def test_label(label):
    """Check label is assigned correctly."""
    qbit = qse.Qbit(label=label)

    if label is None:
        assert qbit.label == "X"
    else:
        assert qbit.label == str(label)


def test_from_qbits():
    """Check Qbit is created correctly from Qbits."""
    nqbit = 4
    qbits = qse.Qbits(positions=np.arange(3 * nqbit).reshape(-1, 3))

    for i in range(nqbit):
        qbit = qse.Qbit(qbits=qbits, index=i)
        assert qbit.label == str(i)
        assert np.allclose(qbit.position, np.array([3 * i, 3 * i + 1, 3 * i + 2]))


def test_set():
    """Check the set method"""
    qbit = qse.Qbit()

    qbit.set("label", "test_qubit")
    assert qbit.label == "test_qubit"

    qbit.set("position", np.ones(3))
    assert np.allclose(qbit.position, np.ones(3))


def test_set_with_qbits():
    """Check the set method words when the Qbit is attached to a Qbits object."""
    nqbit = 4
    qbits = qse.Qbits(positions=np.arange(3 * nqbit).reshape(-1, 3))
    index = 1

    qbit = qse.Qbit(qbits=qbits, index=index)

    label = "test_qubit"
    qbit.set("label", label)
    assert qbit.label == label
    assert qbits.labels[index] == label

    position = 4 * np.ones(3)
    qbit.set("position", position)
    assert np.allclose(qbit.position, position)
    assert np.allclose(qbits.positions[index], position)
