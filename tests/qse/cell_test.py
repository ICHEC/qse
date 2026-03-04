import numpy as np
import pytest

import qse


@pytest.mark.parametrize(
    "cell, expected",
    [
        (4.232 * np.eye(3), np.eye(3) / 4.232),  # square
        (
            np.array(  # rectangular
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 2.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0 / 2.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
        ),
        (
            np.array(  # hex
                [
                    [1.0, 0.0, 0.0],
                    [0.5, 0.5 * np.sqrt(3), 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, -1.0 / np.sqrt(3), 0.0],
                    [0.0, 2.0 / np.sqrt(3), 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
        ),
        (
            0.5
            * np.array(  # fcc
                [
                    [0.0, 1.0, 1.0],
                    [1.0, 0.0, 1.0],
                    [1.0, 1.0, 0.0],
                ]
            ),
            np.array(
                [
                    [-1.0, 1.0, 1.0],
                    [1.0, -1.0, 1.0],
                    [1.0, 1.0, -1.0],
                ]
            ),
        ),
    ],
)
def test_reciprocal_cell(cell, expected):
    """Test the reciprocal method."""
    cell = qse.Cell(cell)
    r_cell = cell.reciprocal() / (2 * np.pi)
    assert np.allclose(r_cell, expected)

    # Check the two cells are orthonormal
    assert np.allclose(cell.lattice_vectors @ r_cell.T, np.eye(3))


@pytest.mark.parametrize("lattice_size", [1.0, 2.2, 3.11, 10.333])
def test_volume(lattice_size):
    """Test the volume method on square lattice."""
    cell = qse.Cell(np.eye(3) * lattice_size)
    assert np.isclose(cell.volume(), lattice_size**3)


@pytest.mark.parametrize(
    "x, y1, y2, z",
    [(2.2, 0.4, 2.1, 1.0), (2.3, -0.4, 2.1, 1.0), (-2.2, -0.4, 2.1, 1.0)],
)
def test_volume_arb(x, y1, y2, z):
    """
    Test the volume method on arbitrary 2d lattice.
    The volume of a cell described by vectors (x, 0, 0),
    (y1, y2, 0) and (0, 0, z) can be easily worked out
    to be x*y2*z.
    """
    cell_np = np.array(
        [
            [x, 0.0, 0.0],
            [y1, y2, 0.0],
            [0, 0.0, z],
        ]
    )
    cell = qse.Cell(cell_np)
    assert np.isclose(cell.volume(), abs(x * y2 * z))
