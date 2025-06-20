import numpy as np
import pytest

import qse


def test_default():
    qbit = qse.Qbit()
    assert qbit.label == "X"
    assert np.allclose(qbit.position, np.zeros(3))
    assert np.allclose(qbit.state, np.array([1, 0]))


@pytest.mark.parametrize("label", ["a" * 3, "123a", 14])
def test_label(label):
    qbit = qse.Qbit(label=label)
    assert qbit.label == str(label)


@pytest.mark.parametrize("state", [3 * np.ones(2), [1.0j, -4]])
def test_state_normalized(state):
    qbit = qse.Qbit(label="123", state=state)
    assert isinstance(qbit, qse.Qbit)
    assert np.isclose(np.linalg.norm(qbit.state), 1.0)


def test_from_qbits():
    nqbit = 4
    qbits = qse.Qbits(positions=np.arange(3*nqbit).reshape(-1, 3))

    for i in range(nqbit):
        qbit = qse.Qbit(qbits=qbits, index=i)
        assert qbit.state == (0, 1)
        assert qbit.label == "Q"
        assert qbit.position == (3*i, 3*i+1, 3*i+2)
