import numpy as np
import pytest

import qse
from qse.lattices import (
    hexagonal_lattice,
    kagome_lattice,
    linear_chain,
    square_lattice,
    triangular_lattice,
)

_lattice_spacings = [0.5, 3.1]
_repeats = [2, 3]


def _lattice_checker(qbits, expected_qbits, lattice_spacing, expected_cellpar):
    assert isinstance(qbits, qse.Qbits)
    assert qbits.nqbits == expected_qbits

    # Check the distance between the first qubit and it's nearest neighbour
    # is the lattice spacing.
    assert np.isclose(qbits.get_all_distances()[0][1:].min(), lattice_spacing)

    # Check the cellpar corresponds to what we expect.
    assert np.allclose(qbits.cell.cellpar(), expected_cellpar)


@pytest.mark.parametrize("lattice_spacing", _lattice_spacings)
@pytest.mark.parametrize("N", _repeats)
def test_linear(lattice_spacing, N):
    qbits = linear_chain(lattice_spacing, N)
    _lattice_checker(qbits, N, lattice_spacing, [lattice_spacing * N, 0, 0, 90, 90, 90])


def test_linear_fail():
    with pytest.raises(
        Exception, match="The repeats must be an integer greater than 1. repeats=1"
    ):
        linear_chain(1.0, 1)


@pytest.mark.parametrize("lattice_spacing", _lattice_spacings)
@pytest.mark.parametrize("N1", _repeats)
@pytest.mark.parametrize("N2", _repeats)
@pytest.mark.parametrize(
    "lattice_func, angles",
    [
        (square_lattice, [90, 90, 90]),
        (triangular_lattice, [90, 90, 60]),
    ],
)
def test_square_and_triangular(lattice_func, lattice_spacing, N1, N2, angles):
    qbits = lattice_func(lattice_spacing, N1, N2)
    _lattice_checker(
        qbits,
        N1 * N2,
        lattice_spacing,
        [lattice_spacing * N1, lattice_spacing * N2, 0] + angles,
    )


@pytest.mark.parametrize("lattice_spacing", _lattice_spacings)
@pytest.mark.parametrize("N1", _repeats)
@pytest.mark.parametrize("N2", _repeats)
def test_hexagonal(lattice_spacing, N1, N2):
    qbits = hexagonal_lattice(lattice_spacing, N1, N2)
    _lattice_checker(
        qbits,
        N1 * N2 * 2,
        lattice_spacing,
        [lattice_spacing * N1 * np.sqrt(3), lattice_spacing * N2 * np.sqrt(3), 0]
        + [90, 90, 60],
    )


@pytest.mark.parametrize("lattice_spacing", _lattice_spacings)
@pytest.mark.parametrize("N1", _repeats)
@pytest.mark.parametrize("N2", _repeats)
def test_kagome(lattice_spacing, N1, N2):
    qbits = kagome_lattice(lattice_spacing, N1, N2)
    _lattice_checker(
        qbits,
        N1 * N2 * 3,
        lattice_spacing,
        [lattice_spacing * N1 * 2, lattice_spacing * N2 * 2, 0] + [90, 90, 60],
    )


@pytest.mark.parametrize(
    "lattice_func",
    [
        square_lattice,
        triangular_lattice,
        hexagonal_lattice,
        kagome_lattice,
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
