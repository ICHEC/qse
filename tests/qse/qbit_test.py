import numpy as np
import pytest

import qse


@pytest.mark.parametrize("label", ["a" * 3, "123a"])
def test_label(label):
    qbit = qse.Qbit(label=label)
    assert qbit.data["label"] == label


@pytest.mark.parametrize("label", [None, 14])
def test_label_default(label):
    qbit = qse.Qbit(label=label)
    assert qbit.data["label"] == "X"


@pytest.mark.parametrize("state", [np.ones(2) * 3, [1.0j, -4]])
def test_state_normalized(state):
    qbit = qse.Qbit(label="123", state=np.array([2, 2.0]))
    assert isinstance(qbit, qse.Qbit)
    assert np.isclose(np.linalg.norm(qbit.data["state"]), 1.0)
