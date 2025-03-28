import numpy as np
import pytest

import qse


@pytest.mark.parametrize("nqbits_1", [1, 2, 3])
@pytest.mark.parametrize("nqbits_2", [1, 2, 3])
def test_add(nqbits_1, nqbits_2):
    """Test adding Qbits together."""
    qbits_1 = qse.Qbits(positions=np.ones((nqbits_1, 3)))
    qbits_2 = qse.Qbits(positions=np.ones((nqbits_2, 3)))
    qbits_12 = qbits_1 + qbits_2

    assert isinstance(qbits_12, qse.Qbits)
    assert qbits_12.nqbits == nqbits_1 + nqbits_2


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
@pytest.mark.parametrize("center_of_mass", [0.0, 1.0])
def test_centre_of_mass(positions, center_of_mass):
    """Test get_center_of_mass, note the positions parametrized all have zero COM."""
    qbits = qse.Qbits(positions=positions + center_of_mass)
    assert np.allclose(qbits.get_center_of_mass(), [center_of_mass] * 3)
