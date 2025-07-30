import matplotlib.pyplot as plt
import numpy as np


def draw(qbits, radius=None, show_labels=False, axes="xy", units=None):
    """
    Visualize the positions of a set of qubits.

    Parameters
    ----------
    qbits: qse.Qbits
        The Qbits object.
    radius: float | str
        A cutoff radius for visualizing bonds.
        Pass 'nearest' to set the radius to the smallest
        distance between the passed qubits.
        If no value is passed the bonds will not be visualized.
    show_labels: bool
        Whether to show the labels of the qubits.
        Defaults to False.
    axes : str
        The 2d axes to visualize the qubits in.
        Must be one of 'xy', 'xz', 'yx', 'yz', 'zx', 'zy'.
        Defaults to 'xy'.
    units : str
        The units of distance.
    """
    if axes not in ['xy', 'xz', 'yx', 'yz', 'zx', 'zy']:
        raise Exception("axes must be one of 'xy', 'xz', 'yx', 'yz', 'zx', 'zy'.")

    x = qbits.positions[:, "xyz".index(axes[0])].copy()
    y = qbits.positions[:, "xyz".index(axes[1])].copy()
    
    draw_bonds = True
    if radius is None:
        draw_bonds = False
    else:
        rij = qbits.get_all_distances()
        min_dist = rij[np.logical_not(np.eye(qbits.nqbits, dtype=bool))].min()
        if radius == "nearest":
            radius = min_dist
        elif min_dist > radius:
            draw_bonds = False

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_aspect("equal")

    if draw_bonds:
        f_tol = 1.01  # fractional tolerance
        nearest_neighbours = (rij <= radius * f_tol)
        np.fill_diagonal(nearest_neighbours, False)
        ii, jj = np.where(nearest_neighbours)
        C = rij[nearest_neighbours]
        C /= C.min()

        ax.quiver(
            x[ii],
            y[ii],
            x[jj] - x[ii],
            y[jj] - y[ii],
            linewidth=1,
            angles="xy",
            scale_units="xy",
            scale=2,
            headaxislength=0,
            headlength=0,
            color="gray",
            alpha=1 / C**3,
        )
    ax.plot(x, y, "o", color="green")

    ax.set_xlabel(axes[0] + f" ({units})" if units is not None else axes[0])
    ax.set_ylabel(axes[1] + f" ({units})" if units is not None else axes[1])

    if show_labels:
        for ind in range(qbits.nqbits):
            ax.text(x[ind], y[ind], s=qbits.labels[ind])


def draw_3d(qbits, radius=None, draw_bonds=True):
    """
    Visualize the positions of a set of qubits.

    Parameters
    ----------
    qbits: qse.Qbits
        The Qbits object.
    radius: float
        A cutoff radius for visualizing bonds.
        Defaults to the smallest distance between the passed qubits.
    draw_bonds: bool
        Whether to show bonds between qubits.
        Defaults to True.
    """
    positions = qbits.positions.copy()
    x, y, z = positions.T

    if draw_bonds:
        rij = qbits.get_all_distances()
        rcut0 = rij[np.logical_not(np.eye(qbits.nqbits, dtype=bool))].min()

        if radius is None:
            rcut = rcut0
        else:
            if rcut0 > radius:
                rcut = radius
                draw_bonds = False
            else:
                rcut = radius

    fig = plt.figure()

    ax = fig.add_subplot(projection="3d")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    if draw_bonds:
        nearest_neighbours = rij <= rcut
        np.fill_diagonal(nearest_neighbours, False)
        ii, jj = np.where(nearest_neighbours)
        X, Y, Z = positions[ii].T
        U, V, W = (positions[jj] - positions[ii]).T
        C = rij[nearest_neighbours]
        C /= C.min()

        ax.quiver(
            X,
            Y,
            Z,
            U,
            V,
            W,
            arrow_length_ratio=0,
            linewidth=0.3,
            color="gray",
            alpha=1 / C**3,
        )
    ax.scatter(x, y, z, "o", color="blue")
