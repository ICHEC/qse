import numpy as np
import pytest

import qse


def _lattice_checker(qbits, expected_qbits, lattice_spacing, expected_cellpar):
    assert isinstance(qbits, qse.Qbits)
    assert qbits.nqbits == expected_qbits

    # Check the distance between the first qubit and it's nearest neighbour
    # is the lattice spacing.
    assert np.isclose(qbits.get_all_distances()[0][1:].min(), lattice_spacing)

    # Check the cellpar corresponds to what we expect.
    assert np.allclose(qbits.cell.cellpar(), expected_cellpar)


@pytest.mark.parametrize("lattice_spacing", [0.5, 3.1])
@pytest.mark.parametrize("N", [2, 3])
def test_linear_lattice(lattice_spacing, N):
    qbits = qse.utils.linear(lattice_spacing, N)
    _lattice_checker(qbits, N, lattice_spacing, [lattice_spacing * N, 0, 0, 90, 90, 90])


@pytest.mark.parametrize("lattice_spacing", [0.5, 3.1])
@pytest.mark.parametrize("N1", [2, 3])
@pytest.mark.parametrize("N2", [2, 3])
@pytest.mark.parametrize(
    "lattice_func, angles",
    [
        (qse.utils.squarelattice, [90, 90, 90]),
        (qse.utils.triangularlattice, [90, 90, 60]),
    ],
)
def test_square_and_triangular_lattices(lattice_func, lattice_spacing, N1, N2, angles):
    qbits = lattice_func(lattice_spacing, N1, N2)
    _lattice_checker(
        qbits,
        N1 * N2,
        lattice_spacing,
        [lattice_spacing * N1, lattice_spacing * N2, 0] + angles,
    )


@pytest.mark.parametrize("lattice_spacing", [0.5, 3.1])
@pytest.mark.parametrize("N1", [2, 3])
@pytest.mark.parametrize("N2", [2, 3])
def test_hexagon_lattice(lattice_spacing, N1, N2):
    qbits = qse.utils.hexagonlattice(lattice_spacing, N1, N2)
    _lattice_checker(
        qbits,
        N1 * N2 * 2,
        lattice_spacing,
        [lattice_spacing * N1 * np.sqrt(3), lattice_spacing * N2 * np.sqrt(3), 0]
        + [90, 90, 60],
    )


@pytest.mark.parametrize("lattice_spacing", [0.5, 3.1])
@pytest.mark.parametrize("N1", [2, 3])
@pytest.mark.parametrize("N2", [2, 3])
def test_kagome_lattice(lattice_spacing, N1, N2):
    qbits = qse.utils.kagomelattice(lattice_spacing, N1, N2)
    _lattice_checker(
        qbits,
        N1 * N2 * 3,
        lattice_spacing,
        [lattice_spacing * N1 * 2, lattice_spacing * N2 * 2, 0] + [90, 90, 60],
    )
