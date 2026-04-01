import numpy as np

import qse


def _check_repeats(repeats, param_name):
    """Validate the repeats parameter."""
    if not isinstance(repeats, int) or repeats < 2:
        raise Exception(
            f"The {param_name} must be an integer greater than 1."
            f" {param_name}={repeats}"
        )


def chain(lattice_spacing: float, repeats: int) -> qse.Qbits:
    """
    Generate a Qbits object in linear chain geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats : int
        The number of repeats.
        Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats, "repeats")
    qbits = qse.Qbits(positions=np.zeros((1, 1)), cell=np.array([[lattice_spacing]]))
    qbits = qbits.repeat((repeats))
    qbits.labels = [str(i) for i in range(qbits.nqbits)]
    return qbits


def square(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits:
    """
    Generate a Qbits object in square lattice geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats_x : int
        The number of repeats in the x direction.
        Must be greater than 1.
    repeats_y : int
        The number of repeats in the y direction.
        Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats_x, "repeats_x")
    _check_repeats(repeats_y, "repeats_y")
    qbits = qse.Qbits(
        positions=np.zeros((1, 2)),
        cell=np.array([[lattice_spacing, 0], [0, lattice_spacing]]),
    )
    qbits = qbits.repeat((repeats_x, repeats_y))
    qbits.labels = [str(i) for i in range(qbits.nqbits)]
    return qbits


def triangular(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits:
    """
    Generate a Qbits object in triangular lattice geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats_x : int
        The number of repeats in the x direction.
        Must be greater than 1.
    repeats_y : int
        The number of repeats in the y direction.
        Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats_x, "repeats_x")
    _check_repeats(repeats_y, "repeats_y")
    cell = np.array([[1, 0], [0.5, np.sqrt(3) / 2]]) * lattice_spacing
    qbits = qse.Qbits(positions=np.zeros((1, 2)), cell=cell)
    qbits = qbits.repeat((repeats_x, repeats_y))
    qbits.labels = [str(i) for i in range(qbits.nqbits)]
    return qbits


def hexagonal(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits:
    """
    Generate a Qbits object in hexagonal lattice geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats_x : int
        The number of repeats in the x direction.
        Must be greater than 1.
    repeats_y : int
        The number of repeats in the y direction.
        Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats_x, "repeats_x")
    _check_repeats(repeats_y, "repeats_y")
    cell = (
        np.array([[3 / 2, np.sqrt(3) / 2], [3 / 2, -np.sqrt(3) / 2]]) * lattice_spacing
    )
    qbits = qse.Qbits(positions=np.array([[0, 0], [lattice_spacing, 0]]), cell=cell)
    qbits = qbits.repeat((repeats_x, repeats_y))
    qbits.labels = [str(i) for i in range(qbits.nqbits)]
    return qbits


def kagome(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits:
    """
    Generate a Qbits object in kagome lattice geometry.

    Parameters
    ----------
    lattice_spacing : float
        The lattice spacing.
    repeats_x : int
        The number of repeats in the x direction.
        Must be greater than 1.
    repeats_y : int
        The number of repeats in the y direction.
        Must be greater than 1.

    Returns
    -------
    Qbits
        The Qbits lattice.
    """
    _check_repeats(repeats_x, "repeats_x")
    _check_repeats(repeats_y, "repeats_y")
    cell = np.array([[2, 0], [1, np.sqrt(3)]]) * lattice_spacing
    qbits = qse.Qbits(
        positions=np.array(
            [
                [0, 0],
                [lattice_spacing, 0],
                [lattice_spacing / 2, np.sqrt(3) * lattice_spacing / 2],
            ]
        ),
        cell=cell,
    )
    qbits = qbits.repeat((repeats_x, repeats_y))
    qbits.labels = [str(i) for i in range(qbits.nqbits)]
    return qbits


def ring(spacing: float, nqbits: int) -> qse.Qbits:
    """
    Generate a Qbits object in a ring geometry.

    Parameters
    ----------
    radius : float
        The spacing between the qubits.
    nqbits : int
        Number of qubits on the ring.

    Returns
    -------
    Qbits
        The Qbits object.
    """
    # Using cosine rule
    # spacing^2 = 2radius^2 [ 1 - cos(2pi/nqubits) ]
    radius = spacing * np.sqrt(0.5 / (1.0 - np.cos(2 * np.pi / nqbits)))
    theta = np.arange(nqbits, dtype=float)
    theta *= 2.0 * np.pi / nqbits
    positions = radius * np.array([np.cos(theta), np.sin(theta)]).T
    return qse.Qbits(positions=positions)


def torus(
    n_outer: int, n_inner: int, inner_radius: float, outer_radius: float
) -> qse.Qbits:
    """
    Generate a Qbits object in a torus geometry.

    Parameters
    ----------
    n_outer : int
        Number of points in larger (outer) dimension.
    n_inner : int
        Number of points in smaller (inner) dimension.
    inner_radius : float
        The inner radius.
    outer_radius : float
        The outer radius.

    Returns
    -------
    Qbits
        The Qbits object.
    """
    theta = (2.0 * np.pi / n_outer) * np.arange(n_outer)
    phi = (2.0 * np.pi / n_inner) * np.arange(n_inner)

    positions = (
        np.array(
            [
                np.outer(outer_radius + inner_radius * np.cos(theta), np.cos(phi)),
                np.outer(outer_radius + inner_radius * np.cos(theta), np.sin(phi)),
                np.outer(inner_radius * np.sin(theta), np.ones(phi.shape)),
            ]
        )
        .reshape(3, n_outer * n_inner)
        .T
    )
    return qse.Qbits(positions=positions)
