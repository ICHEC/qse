import numpy as np
import pytest

import qse


def test_empty_qbits():
    qbits = qse.Qbits()
    assert isinstance(qbits, qse.Qbits)


@pytest.mark.parametrize("labels", [None, ["a", "b", "c"], ["qb1", "interesting_qbit"]])
def test_labels(labels):
    qbits = qse.Qbits(labels=labels)
    assert isinstance(qbits, qse.Qbits)


@pytest.mark.parametrize("labels", [4, "hello", "x"])
def test_labels_fail(labels):
    with pytest.raises(Exception, match="'labels' must be a list."):
        qse.Qbits(labels=labels)


@pytest.mark.parametrize("positions", [np.arange(9).reshape(-1, 3), np.zeros((1, 3))])
def test_positions(positions):
    qbits = qse.Qbits(positions=positions)
    assert isinstance(qbits, qse.Qbits)


@pytest.mark.parametrize("positions", [np.zeros(3), [1, 2]])
def test_positions_fail(positions):
    with pytest.raises(Exception):
        qse.Qbits(positions=positions)


def test_positions_and_scaled_fail():
    with pytest.raises(
        Exception,
        match="Both 'positions' and 'scaled_positions' cannot be passed at the same time.",
    ):
        qse.Qbits(
            positions=np.arange(9).reshape(-1, 3),
            scaled_positions=np.arange(9).reshape(-1, 3),
        )


def test_scaled_positions_fail():
    with pytest.raises(
        Exception,
        match="'scaled_positions' requires 'cell' to not be None.",
    ):
        qse.Qbits(
            scaled_positions=np.arange(9).reshape(-1, 3),
        )


def test_labels_positions_fail():
    positions = np.arange(9).reshape(-1, 3)
    labels = ["a", "b"]
    with pytest.raises(
        Exception, match="Both 'positions' and 'labels' must have the same length."
    ):
        qse.Qbits(positions=positions, labels=labels)


def test_init_from_qbit():
    position = [0, 0, 0]
    label = "test_qbit"
    qbit = qse.Qbit(position=position, label=label)
    qbits = qse.Qbits.from_qbit_list([qbit])
    assert isinstance(qbits, qse.Qbits)
    assert qbits.nqbits == 1
    assert qbits.labels == [label]
    assert np.allclose(qbits.positions, [position])
