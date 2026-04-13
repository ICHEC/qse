"""
Visualize
---------

This is the visualisation module of QSE package.
It defines a number of helper functions to better
visualise data.

One of the most common function is the `draw` function
to draw the qubits.

"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

rads = np.linspace(12, 1, 10)
colors = np.sqrt(np.linspace(0.1, 0.9, 10))


def draw(
    qbits, radius=None, show_labels=False, colouring=None, units=None, equal_aspect=True
):
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
    colouring: str | list
        A set of integers used to assign different colors to each Qubit.
        This can be used to view different magnetic orderings.
        Must have the same length as the number of Qubits.
    units : str, optional
        The units of distance.
    equal_aspect : bool, optional
        Whether to have the same scaling for the axes.
        Defaults to True.
    """
    if colouring is not None:
        if len(colouring) != qbits.nqbits:
            raise Exception("The length of colouring must equal the number of Qubits.")
        colouring = [int(i) for i in colouring]

    draw_bonds = False if radius is None else True
    rij = None
    min_dist = None

    if draw_bonds:
        rij = qbits.get_all_distances()
        min_dist = rij[np.logical_not(np.eye(qbits.nqbits, dtype=bool))].min()
        if radius == "nearest":
            radius = min_dist
        elif min_dist > radius:
            draw_bonds = False

    fig = plt.figure()
    projection = "3d" if qbits.dim == 3 else None
    ax = fig.add_subplot(projection=projection)
    if equal_aspect:
        ax.set_aspect("equal")

    if qbits.dim == 3:
        _draw_3d(qbits, draw_bonds, radius, rij, min_dist, ax)
    else:
        _draw_2d(
            qbits,
            draw_bonds,
            radius,
            rij,
            min_dist,
            units,
            colouring,
            show_labels,
            ax,
        )
    return fig


def _draw_3d(qbits, draw_bonds, radius, rij, min_dist, ax):
    positions = qbits.positions

    if draw_bonds:
        f_tol = 1.01  # fractional tolerance
        neighbours = rij <= radius * f_tol
        np.fill_diagonal(neighbours, False)
        ii, jj = np.where(neighbours)
        X, Y, Z = positions[ii].T
        U, V, W = (positions[jj] - positions[ii]).T
        alpha = (min_dist / rij[neighbours]) ** 3

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
            alpha=alpha,
        )
    x, y, z = positions.T

    for r, c in zip(rads, colors):
        ax.scatter(x, y, z, s=r**2, color=(0.1, c, 0.5), zorder=1, alpha=0.8)


def _draw_2d(
    qbits,
    draw_bonds,
    radius,
    rij,
    min_dist,
    units,
    colouring,
    show_labels,
    ax,
):
    ax.set_xlabel("x" + f" ({units})" if units is not None else "x")
    ax.set_ylabel("y" + f" ({units})" if units is not None else "y")

    if qbits.dim == 2:
        x, y = qbits.positions.T
    else:
        x = qbits.positions.T.flatten()
        y = np.zeros(qbits.nqbits)

    if draw_bonds:
        f_tol = 1.01  # fractional tolerance
        neighbours = np.array(
            [
                (i, j)
                for i in range(qbits.nqbits - 1)
                for j in range(i + 1, qbits.nqbits)
                if rij[i, j] <= radius * f_tol
            ]
        )
        for i, j in neighbours:
            alpha = (min_dist / rij[i, j]) ** 3
            ax.plot([x[i], x[j]], [y[i], y[j]], c="gray", alpha=alpha, zorder=-1)

    if colouring is not None:
        inds0 = [j == 0 for j in colouring]
        inds1 = [j == 1 for j in colouring]

        for r, c in zip(rads, colors):
            ax.scatter(
                x[inds0], y[inds0], s=r**2, color=(0.1, c, 0.5), zorder=1, alpha=0.8
            )
            ax.scatter(
                x[inds1], y[inds1], s=r**2, color=(c, 0.1, 0.5), zorder=1, alpha=0.8
            )

    else:
        for r, c in zip(rads, colors):
            ax.scatter(x, y, s=r**2, color=(0.1, c, 0.5), zorder=1, alpha=0.8)

    if show_labels:
        for ind in range(qbits.nqbits):
            ax.text(x[ind], y[ind], s=qbits.labels[ind])


def view_matrix(matrix, labels_x=None, labels_y=None, vcenter=None):
    """
    Visualise a matrix.

    Parameters
    ----------
    matrix: np.ndarray
        The matrix to be visualised.
    labels_x: list, optional
        Labels to be displayed on the x axis.
    labels_y: list, optional
        Labels to be displayed on the y axis.
    vcenter: float, optional
        The center of the colorbar.
    """

    fig = plt.figure()
    ax = fig.add_subplot()

    if labels_x is not None:
        ax.set_xticks(range(matrix.shape[0]), labels_x)

    if labels_y is not None:
        ax.set_yticks(range(matrix.shape[1]), labels_y)

    norm = None
    if vcenter is not None:
        norm = TwoSlopeNorm(vmin=matrix.min(), vcenter=vcenter, vmax=matrix.max())

    im = ax.imshow(
        matrix,
        cmap="PiYG",
        norm=norm,
    )

    fig.colorbar(im, ax=ax)

    ax.invert_yaxis()  # More natural to invert y-axis for these plots.
    return fig


def bar(dict, cutoff=0, ylabel="Count"):
    """
    Plot a bar chart from a dictionary, filtering values below a cutoff.

    Parameters
    ----------
    dict : dict
        A dictionary where keys are categories and values are their corresponding
        counts or values.
    cutoff : int or float, optional
        Minimum value threshold for inclusion in the plot. Keys with values less
        than or equal to `cutoff` are excluded.
        Default is 0.
    ylabel : str, optional
        Label for the y-axis of the plot. Default is "Count".
    """
    fig = plt.figure()

    dict = {k: v for k, v in dict.items() if v > cutoff}
    plt.bar(list(dict.keys()), list(dict.values()))
    plt.xticks(rotation="vertical")
    plt.ylabel(ylabel)

    return fig
