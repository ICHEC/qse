import numpy as np
import pytest

import qse
from qse.lattices import (
    chain,
    hexagonal,
    kagome,
    ring,
    square,
    triangular,
)

_lattice_spacings = [0.5, 3.1]
_repeats = [2, 3]


def _lattice_checker(
    qbits, expected_qbits, lattice_spacing, expected_cell, expected_positions
):
    assert isinstance(qbits, qse.Qbits)
    assert qbits.nqbits == expected_qbits
    assert qbits.positions.shape == (expected_qbits, 3)

    # Check the labels
    assert all(qbits.labels[i] == str(i) for i in range(qbits.nqbits))

    # Check the distance between the first qubit and it's nearest neighbour
    # is the lattice spacing.
    assert np.isclose(qbits.get_all_distances()[0][1:].min(), lattice_spacing)

    # Check the cell
    assert np.allclose(qbits.cell, expected_cell)

    # Check the positions
    assert np.allclose(qbits.positions, expected_positions)


@pytest.mark.parametrize("lattice_spacing", _lattice_spacings)
@pytest.mark.parametrize("N", _repeats)
def test_linear(lattice_spacing, N):
    qbits = chain(lattice_spacing, N)
    expected_cell = np.zeros((3, 3))
    expected_cell[0, 0] = lattice_spacing * N
    expected_positions = lattice_spacing * np.array([[i, 0, 0] for i in range(N)])
    _lattice_checker(qbits, N, lattice_spacing, expected_cell, expected_positions)


def test_linear_fail():
    with pytest.raises(
        Exception, match="The repeats must be an integer greater than 1. repeats=1"
    ):
        chain(1.0, 1)


@pytest.mark.parametrize("lattice_spacing", _lattice_spacings)
@pytest.mark.parametrize("N1", _repeats)
@pytest.mark.parametrize("N2", _repeats)
def test_square(lattice_spacing, N1, N2):
    qbits = square(lattice_spacing, N1, N2)
    expected_cell = np.zeros((3, 3))
    expected_cell[0, 0] = lattice_spacing * N1
    expected_cell[1, 1] = lattice_spacing * N2
    expected_positions = lattice_spacing * np.array(
        [[i, j, 0.0] for i in range(N1) for j in range(N2)]
    )
    _lattice_checker(
        qbits,
        N1 * N2,
        lattice_spacing,
        expected_cell,
        expected_positions,
    )


@pytest.mark.parametrize("lattice_spacing", _lattice_spacings)
@pytest.mark.parametrize("N1", _repeats)
@pytest.mark.parametrize("N2", _repeats)
def test_triangular(lattice_spacing, N1, N2):
    qbits = triangular(lattice_spacing, N1, N2)

    expected_cell = lattice_spacing * np.array(
        [
            [N1, 0.0, 0.0],
            [N2 * np.cos(np.pi / 3), N2 * np.sin(np.pi / 3), 0.0],
            [0.0, 0.0, 0.0],
        ]
    )

    expected_positions = lattice_spacing * np.array(
        [
            [(i + j * 0.5), j * np.sqrt(3) * 0.5, 0.0]
            for i in range(N1)
            for j in range(N2)
        ]
    )
    _lattice_checker(
        qbits,
        N1 * N2,
        lattice_spacing,
        expected_cell,
        expected_positions,
    )


@pytest.mark.parametrize("lattice_spacing", _lattice_spacings)
@pytest.mark.parametrize("N1", _repeats)
@pytest.mark.parametrize("N2", _repeats)
def test_hexagonal(lattice_spacing, N1, N2):
    qbits = hexagonal(lattice_spacing, N1, N2)

    expected_cell = (
        lattice_spacing
        * np.sqrt(3)
        * 0.5
        * np.array(
            [
                [N1 * np.sqrt(3), N1, 0.0],
                [N2 * np.sqrt(3), -N2, 0.0],
                [0.0, 0.0, 0.0],
            ]
        )
    )

    arr1 = lattice_spacing * np.array(
        [
            [(i + j) * 1.5, (i - j) * np.sqrt(3) * 0.5, 0.0]
            for i in range(N1)
            for j in range(N2)
        ]
    )
    arr2 = arr1.copy()
    arr2[:, 0] += lattice_spacing
    expected_positions = np.zeros((len(arr1) + len(arr2), 3))
    expected_positions[0::2] = arr1
    expected_positions[1::2] = arr2
    _lattice_checker(
        qbits,
        N1 * N2 * 2,
        lattice_spacing,
        expected_cell,
        expected_positions,
    )


@pytest.mark.parametrize("lattice_spacing", _lattice_spacings)
@pytest.mark.parametrize("N1", _repeats)
@pytest.mark.parametrize("N2", _repeats)
def test_kagome(lattice_spacing, N1, N2):
    qbits = kagome(lattice_spacing, N1, N2)

    expected_cell = lattice_spacing * np.array(
        [
            [N1 * 2, 0.0, 0.0],
            [N2, N2 * np.sqrt(3), 0.0],
            [0.0, 0.0, 0.0],
        ]
    )

    arr1 = lattice_spacing * np.array(
        [[(2 * i + j), j * np.sqrt(3), 0.0] for i in range(N1) for j in range(N2)]
    )
    arr2 = arr1.copy()
    arr2[:, 0] += lattice_spacing
    arr3 = arr1.copy()
    arr3[:, 0] += lattice_spacing * 0.5
    arr3[:, 1] += np.sqrt(3) * lattice_spacing / 2
    expected_positions = np.zeros((len(arr1) + len(arr2) + len(arr3), 3))
    expected_positions[0::3] = arr1
    expected_positions[1::3] = arr2
    expected_positions[2::3] = arr3

    _lattice_checker(
        qbits,
        N1 * N2 * 3,
        lattice_spacing,
        expected_cell,
        expected_positions,
    )


@pytest.mark.parametrize(
    "lattice_func",
    [
        square,
        triangular,
        hexagonal,
        kagome,
    ],
)
def test_lattices_fail(lattice_func):
    with pytest.raises(
        Exception, match="The repeats_x must be an integer greater than 1. repeats_x=1"
    ):
        lattice_func(1.0, 1, 2)

    with pytest.raises(
        Exception, match="The repeats_y must be an integer greater than 1. repeats_y=1"
    ):
        lattice_func(1.0, 2, 1)


@pytest.mark.parametrize("spacing", [0.5, 2.98])
@pytest.mark.parametrize("nqbits", [2, 4, 7])
def test_ring(spacing, nqbits):
    """Test the nearest neighbour distance is the spacing in the ring geometry."""
    qbits = ring(spacing, nqbits)
    distances = qbits.get_all_distances()
    assert np.allclose(np.diag(distances, k=1), spacing)
    assert np.allclose(np.diag(distances, k=-1), spacing)
