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


@pytest.mark.parametrize("state", [np.ones(2) * 3, [1.0j, -4]])
def test_state_normalized(state):
    qbit = qse.Qbit(label="123", state=state)
    assert isinstance(qbit, qse.Qbit)
    assert np.isclose(np.linalg.norm(qbit.state), 1.0)
