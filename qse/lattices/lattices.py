import numpy as np

import qse


def _check_repeats(repeats, param_name):
    """Validate the repeats parameter."""
    if not isinstance(repeats, int) or repeats < 2:
        raise Exception(
            f"The {param_name} must be an integer greater than 1."
            f" {param_name}={repeats}"
        )


def _lattice_creator(
    unit_cell_2d: list,
    lattice_spacing: float,
    repeats: tuple[int],
    qubit_positions: list = [[0, 0, 0]],
):
    """
    Create a Qbits object from a 2d cell, lattice_spacing,
    repeats and qubit positions.
    """
    unit = qse.Qbits(positions=np.array(qubit_positions))

    # We add [[0, 0, 0]] to convert the unit cell into 3d.
    unit.cell = lattice_spacing * qse.cell.Cell(unit_cell_2d + [[0, 0, 0]])
    return unit.repeat(repeats + (1,))


def chain(lattice_spacing: float = 1.0, repeats: int = 6):
    """
    Generate a Qbits object in linear chain geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats : int
        The number of repeats. Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats, "repeats")
    return _lattice_creator([[1, 0, 0], [0, 0, 0]], lattice_spacing, (repeats, 1))


def square(lattice_spacing: float = 1.0, repeats_x: int = 2, repeats_y: int = 2):
    """
    Generate a Qbits object in square lattice geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats_x : int
        The number of repeats in the x direction. Must be greater than 1.
    repeats_y : int
        The number of repeats in the y direction. Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats_x, "repeats_x")
    _check_repeats(repeats_y, "repeats_y")
    return _lattice_creator(
        [[1, 0, 0], [0, 1, 0]], lattice_spacing, (repeats_x, repeats_y)
    )


def triangular(lattice_spacing: float = 1.0, repeats_x: int = 2, repeats_y: int = 2):
    """
    Generate a Qbits object in triangular lattice geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats_x : int
        The number of repeats in the x direction. Must be greater than 1.
    repeats_y : int
        The number of repeats in the y direction. Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats_x, "repeats_x")
    _check_repeats(repeats_y, "repeats_y")
    return _lattice_creator(
        [[1, 0, 0], [0.5, np.sqrt(3) / 2, 0]], lattice_spacing, (repeats_x, repeats_y)
    )


def hexagonal(lattice_spacing: float = 1.0, repeats_x: int = 2, repeats_y: int = 2):
    """
    Generate a Qbits object in hexagonal lattice geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats_x : int
        The number of repeats in the x direction. Must be greater than 1.
    repeats_y : int
        The number of repeats in the y direction. Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats_x, "repeats_x")
    _check_repeats(repeats_y, "repeats_y")
    return _lattice_creator(
        [
            [3 / 2, np.sqrt(3) / 2, 0],
            [3 / 2, -np.sqrt(3) / 2, 0],
        ],
        lattice_spacing,
        (repeats_x, repeats_y),
        [[0, 0, 0], [lattice_spacing, 0, 0]],
    )


def kagome(lattice_spacing: float = 1.0, repeats_x: int = 2, repeats_y: int = 2):
    """
    Generate a Qbits object in kagome lattice geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats_x : int
        The number of repeats in the x direction. Must be greater than 1.
    repeats_y : int
        The number of repeats in the y direction. Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats_x, "repeats_x")
    _check_repeats(repeats_y, "repeats_y")
    return _lattice_creator(
        [[2, 0, 0], [1, np.sqrt(3), 0]],
        lattice_spacing,
        (repeats_x, repeats_y),
        [
            [0, 0, 0],
            [lattice_spacing, 0, 0],
            [lattice_spacing / 2, np.sqrt(3) * lattice_spacing / 2, 0],
        ],
    )
