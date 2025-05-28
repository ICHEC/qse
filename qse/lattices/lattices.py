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
    """Create a Qbits object from a 2d cell,
    lattice_spacing, repeats and qubit positions.

    Parameters
    ----------
    unit_cell_2d : list
        The unit cell
    lattice_spacing : float
        Lattice spacing
    repeats : tuple[int]
        Repetitions along each direction, such as (3, 4, 2)
    qubit_positions : list, optional
        Positions for qubits as Nx3 array, by default [[0, 0, 0]]

    Returns
    -------
    Qbits
        Qbits object repeated along desired directions.
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


def ring(center=np.zeros(3), radius=3.0, npoints=12):
    """Generate Qbits object in ring geometry.
        The ring is placed in xy plane at center.

    Parameters
    ----------
    center : array_like, optional
        3D coordinate for the centre
        of the ring, by default np.zeros(3)
    radius : float, optional
        Radius of ring, by default 3.0
    npoints : int, optional
        Number of points on ring, by default 12

    Returns
    -------
    Qbits
        The Qbits object with ring geometry.
    """
    theta = np.arange(npoints, dtype=float)
    theta *= 2.0 * np.pi / npoints
    positions = np.array([np.cos(theta), np.sin(theta), np.zeros(npoints)]).T
    positions *= radius
    positions += center
    qring = qse.Qbits(positions=positions)
    return qring


def torus(N1=12, N2=12, Rin=1.0, Rout=3.0, center=np.zeros(3)):
    """Generate Qbits object in a torus geometry.

    Parameters
    ----------
    N1 : int, optional
        Number of points in larger (outer) dimension, by default 12
    N2 : int, optional
        Number of points in smaller (inner) dimension, by default 12
    Rin : float, optional
        Inner radius, by default 1.0
    Rout : float, optional
        Outer radius, by default 3.0
    center : _type_, optional
        Center of the torus, by default np.zeros(3)

    Returns
    -------
    Qbits
        The Qbits object with torus geometry.
    """
    theta = (2.0 * np.pi / N1) * np.arange(N1)
    phi = (2.0 * np.pi / N2) * np.arange(N2)
    #
    positions = (
        np.array(
            [
                np.outer(Rout + Rin * np.cos(theta), np.cos(phi)),
                np.outer(Rout + Rin * np.cos(theta), np.sin(phi)),
                np.outer(Rin * np.sin(theta), np.ones(phi.shape)),
            ]
        )
        .reshape(3, N1 * N2)
        .T
    )
    positions += center
    qtorus = qse.Qbits(positions=positions)
    return qtorus
