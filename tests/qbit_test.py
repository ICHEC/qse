import numpy as np
import pytest

import qse


def test_default():
    """Check default behaviour."""
    qbit = qse.Qbit()
    assert isinstance(qbit, qse.Qbit)
    assert qbit.label == "X"
    assert np.allclose(qbit.position, np.zeros(3))
    assert np.allclose(qbit.state, np.array([1, 0]))


@pytest.mark.parametrize("label", ["a" * 3, "123a", None, 14])
def test_label(label):
    """Check label is assigned correctly."""
    qbit = qse.Qbit(label=label)
    if not isinstance(label, str):
        # default behaviour
        assert qbit.label == "X"
    else:
        assert qbit.label == label


@pytest.mark.parametrize("state", [np.ones(2) * 3, [1.0j, -4]])
def test_state(state):
    """Check state is normalized."""
    qbit = qse.Qbit(label="123", state=state)
    assert isinstance(qbit, qse.Qbit)
    assert np.isclose(np.linalg.norm(qbit.state), 1.0)


def test_from_qbits():
    """Check Qbit is created correctly from Qbits."""
    nqbit = 4
    qbits = qse.Qbits(positions=np.arange(3 * nqbit).reshape(-1, 3))

    for i in range(nqbit):
        qbit = qse.Qbit(qbits=qbits, index=i)
        assert np.allclose(qbit.state, np.array([1, 0]))
        assert qbit.label == str(i)
        assert np.allclose(qbit.position, np.array([3 * i, 3 * i + 1, 3 * i + 2]))


def test_set():
    """Check the set method"""
    qbit = qse.Qbit()
    qbit.set("state", np.array([1.0, 1.0j]) / np.sqrt(2))
    assert np.allclose(qbit.state, np.array([1.0, 1.0j]) / np.sqrt(2))

    qbit.set("label", "test_qubit")
    assert qbit.label == "test_qubit"

    qbit.set("position", np.ones(3))
    assert np.allclose(qbit.position, np.ones(3))


def test_from_set_qbits():
    """Check the set method words when the Qbit is attached to a Qbits object."""
    nqbit = 4
    qbits = qse.Qbits(positions=np.arange(3 * nqbit).reshape(-1, 3))
    index = 1

    qbit = qse.Qbit(qbits=qbits, index=index)
    qbit.set("state", np.array([1.0, 1.0j]) / np.sqrt(2))
    assert np.allclose(qbit.state, np.array([1.0, 1.0j]) / np.sqrt(2))
    assert np.allclose(qbits.states[index], np.array([1.0, 1.0j]) / np.sqrt(2))

    qbit.set("label", "g")
    assert qbit.label == "g"
    assert qbits.labels[index] == "g"

    qbit.set("position", np.ones(3))
    assert np.allclose(qbit.position, np.ones(3))
    assert np.allclose(qbits.positions[index], np.ones(3))
