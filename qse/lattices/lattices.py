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
    Generate a Qbits object in linear chain geometry. This is
    the simplest lattice in one dimension, with equally separated
    qubits along a line.

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

    Examples
    ----
    .. jupyter-execute::

        import qse
        q1d = qse.lattices.chain(
            lattice_spacing=2.0,
            repeats=6)
        q1d.draw()

    Note
    ----
    For a one dimensional lattice, there isn't really a 'lattice vector'
    as the positions are represented as numbers. The lattice separation
    :math:`a` can be thought of as 1D lattice vector.
    """
    _check_repeats(repeats, "repeats")
    qbits = qse.Qbits(positions=np.zeros((1, 1)), cell=np.array([[lattice_spacing]]))
    qbits = qbits.repeat((repeats))
    qbits.labels = [str(i) for i in range(qbits.nqbits)]
    return qbits


def square(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits:
    r"""
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


    Examples
    ----
    .. jupyter-execute::

        import qse
        qbit = qse.lattices.square(
            lattice_spacing=2.0,
            repeats_x=3, repeats_y=3)
        qbit.draw(radius=4)

    Note
    ----
    For the square lattice the unit vectors are

    .. math::
        A_1 = (a, 0) \quad A_2 = (0, a)

    with the single basis point at origin :math:`(0, 0)`.
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
    r"""
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

    Examples
    ----
    .. jupyter-execute::

        import qse
        qbit = qse.lattices.triangular(
            lattice_spacing=2.0,
            repeats_x=3, repeats_y=3)
        qbit.draw(radius=3)

    Note
    -----
    For the triangular lattice the unit vectors are

    .. math::
        A_1 = (a, 0) \quad A_2 = (\frac 12 a,\frac{\sqrt{3}}{2} a)

    with the single basis point at origin :math:`(0, 0)`.
    """
    _check_repeats(repeats_x, "repeats_x")
    _check_repeats(repeats_y, "repeats_y")
    cell = np.array([[1, 0], [0.5, np.sqrt(3) / 2]]) * lattice_spacing
    qbits = qse.Qbits(positions=np.zeros((1, 2)), cell=cell)
    qbits = qbits.repeat((repeats_x, repeats_y))
    qbits.labels = [str(i) for i in range(qbits.nqbits)]
    return qbits


def hexagonal(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits:
    r"""
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

    Examples
    ----
    .. jupyter-execute::

        import qse
        qbit = qse.lattices.hexagonal(
            lattice_spacing=2.0,
            repeats_x=3, repeats_y=3)
        qbit.draw(radius=3)

    Note
    -----
    For the hexagonal lattice the unit vectors are

    .. math::
        A_1 = (\frac 32 a, \frac{\sqrt{3}}{2} a)\quad
        A_2 = (\frac 32 a, -\frac{\sqrt{3}}{2} a)

    with two basis points at :math:`(0, 0)` and :math:`(a, 0)`.
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
    r"""
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

    Examples
    ----
    .. jupyter-execute::

        import qse
        qbit = qse.lattices.kagome(
            lattice_spacing=2.0,
            repeats_x=3, repeats_y=3)
        qbit.draw(radius=3)

    Note
    -----
    For the kagome lattice the unit vectors are

    .. math::
        A_1 = (2a, 0)\quad
        A_2 = (a, \sqrt{3}a)

    with three basis points at

    .. math::
        (0, 0)\quad (a, 0) \quad (\frac 12 a, \frac{\sqrt{3}}{2}a)
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
    r"""
    Generate a Qbits object in a ring geometry.

    Parameters
    ----------
    spacing : float
        The spacing between the qubits.
    nqbits : int
        Number of qubits on the ring.

    Returns
    -------
    Qbits
        The Qbits object.

    Examples
    --------
    .. jupyter-execute::

        import qse
        qbit = qse.lattices.ring(
            spacing=1.0,
            nqbits=9)
        qbit.draw(radius=3)

    Note
    ----
    This can be thought of as 1D chain where you connect
    the end point and place them on a circular shape.

    :math:`n` points placed on a circle of radius :math:`r`
    have coordinates :math:`X_l = (r\cos\theta_l, r\sin\theta_l)`,
    where :math:`l` ranges from 0 to :math:`n-1` and
    :math:`\theta_l = \frac{2\pi l}{n}`.

    If the qubit spacing supplied is :math:`a`, then we compute the
    appropriate radius by solving :math:`a = |X_{l+1} - X_l|`, which
    gives us the radius :math:`r` as

    .. math::
        r = \frac{a} {\sqrt{2 (1 - \cos{\frac{2\pi}{n}})}}
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

    Examples
    --------
    .. jupyter-execute::

        import qse
        qbit = qse.lattices.torus(
            n_outer=12, n_inner=12,
            inner_radius=5.0,
            outer_radius=6.0)
        qbit.draw()

    Note
    ----
    Though its not a commonly explored geometry, it represents
    the 2D version of a ring, with periodic connection built in.
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
