import matplotlib.pyplot as plt
import numpy as np

rads = np.linspace(12, 1, 10)
colors = np.sqrt(np.linspace(0.1, 0.9, 10))


def draw(qbits, radius=None, show_labels=False, colouring=None, units=None):
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
    """
    if colouring is not None:
        if len(colouring) != qbits.nqbits:
            raise Exception("The length of colouring must equal the number of Qubits.")
        colouring = [int(i) for i in colouring]

    cell_rank = qbits.cell.rank
    position_rank = np.linalg.matrix_rank(qbits.positions)
    # if cell_rank is more than position_rank, it means that a higher dimensional cell
    # is present, and actual structure requires repettition of cells to
    # clearly visualize. if position_rank is more than cell_rank, it means the repettion
    # happens in lower dimension of a local geometry which visualized in
    # higher dimension. Either way, we need to see things in higher dimension.

    rank = max(cell_rank, position_rank)

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

    if rank == 3:
        _draw_3d(qbits, draw_bonds, radius, rij, min_dist)
    else:
        _draw_2d(
            qbits, draw_bonds, radius, rij, min_dist, units, colouring, show_labels
        )


def _draw_3d(qbits, draw_bonds, radius, rij, min_dist):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.set_aspect("equal")

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

    return fig


def _draw_2d(qbits, draw_bonds, radius, rij, min_dist, units, colouring, show_labels):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_aspect("equal")
    ax.set_xlabel("x" + f" ({units})" if units is not None else "x")
    ax.set_ylabel("y" + f" ({units})" if units is not None else "y")

    x, y, _ = qbits.positions.T
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

    return fig
