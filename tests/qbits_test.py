import numpy as np
import pytest

import qse


def _qbits_checker(qbits, positions, total_qubits):
    assert isinstance(qbits, qse.Qbits)
    assert np.allclose(qbits.get_positions(), positions)
    assert qbits.nqbits == total_qubits


@pytest.mark.parametrize("nqbits", [1, 2, 3, 10])
def test_basic_properties(nqbits):
    """Test basic Qbits methods."""
    positions = np.random.rand(nqbits, 3)
    qbits = qse.Qbits(positions=positions)
    _qbits_checker(qbits, positions, nqbits)
    assert len(qbits) == nqbits


@pytest.mark.parametrize("nqbits_1", [1, 2, 3])
@pytest.mark.parametrize("nqbits_2", [1, 2, 3])
def test_add(nqbits_1, nqbits_2):
    """Test adding Qbits together."""
    positions_1 = np.random.rand(nqbits_1, 3)
    qbits_1 = qse.Qbits(positions=positions_1)

    positions_2 = np.random.rand(nqbits_2, 3)
    qbits_2 = qse.Qbits(positions=positions_2)

    qbits_12 = qbits_1 + qbits_2

    _qbits_checker(
        qbits_12, np.concatenate((positions_1, positions_2)), nqbits_1 + nqbits_2
    )


def test_get_all_distances():
    """Test get_all_distances on a simple set of qbits."""
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
    )

    qbits = qse.Qbits(positions=positions)

    assert qbits.nqbits == positions.shape[0]

    distances = np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, np.sqrt(2)],
            [1.0, np.sqrt(2), 0.0],
        ]
    )

    assert np.allclose(qbits.get_all_distances(), distances)


@pytest.mark.parametrize(
    "positions",
    [
        np.array([[-1, 0, 0], [1, 0, 0]]),
        np.array([[-1, 0, 0], [1, 0, 0], [1, 2, 0], [1, -2, 6], [-2, 0, -6]]),
    ],
)
@pytest.mark.parametrize("centroid", [0.0, 1.0, -1.0])
def test_get_centroid(positions, centroid):
    """Test get_centroid, note the positions parametrized all have zero centroid."""
    qbits = qse.Qbits(positions=positions + centroid)
    assert np.allclose(qbits.get_centroid(), [centroid] * 3)


@pytest.mark.parametrize(
    "positions",
    [
        np.array([[-1, 0, 0], [1, 0, 0]]),
        np.array([[-1, 0, 0], [1, 0, 0], [1, 2, 0], [1, -2, 6], [-2, 0, -6]]),
    ],
)
@pytest.mark.parametrize(
    "centroid", [[0.0] * 3, [1.0] * 3, [-1.0] * 3, np.random.rand(3)]
)
def test_set_centroid(positions, centroid):
    """Test set_centroid, note the positions parametrized all have zero centroid."""
    qbits = qse.Qbits(positions=positions)
    qbits.set_centroid(centroid)
    assert np.allclose(qbits.get_centroid(), centroid)


@pytest.mark.parametrize("nqbits", [1, 2, 3, 10])
def test_rattle(nqbits):
    """Test rattle."""
    positions = np.random.rand(nqbits, 3)
    qbits = qse.Qbits(positions=positions)
    qbits.rattle()

    assert qbits.get_positions().shape == positions.shape
    assert not np.allclose(qbits.get_positions(), positions)
