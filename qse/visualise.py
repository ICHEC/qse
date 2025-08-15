import matplotlib.pyplot as plt
import numpy as np


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
    if (colouring is not None) and (len(colouring) != qbits.nqbits):
        raise Exception("The length of colouring must equal the number of Qubits.")

    cell_rank = qbits.cell.rank
    position_rank = np.linalg.matrix_rank(qbits.positions)
    # if cell_rank is more than position_rank, it means that a higher dimensional cell
    # is present, and actual structure requires repettition of cells to
    # clearly visualize. if position_rank is more than cell_rank, it means the repettion
    # happens in lower dimension of a local geometry which visualized in
    # higher dimension. Either way, we need to see things in higher dimension.

    rank = max(cell_rank, position_rank)
    positions = qbits.positions.copy()
    x, y, z = positions.T

    draw_bonds = False if radius is None else True

    if draw_bonds:
        rij = qbits.get_all_distances()
        min_dist = rij[np.logical_not(np.eye(qbits.nqbits, dtype=bool))].min()
        if radius == "nearest":
            radius = min_dist
        elif min_dist > radius:
            draw_bonds = False

    if draw_bonds:
        f_tol = 1.01  # fractional tolerance
        neighbours = rij <= radius * f_tol
        np.fill_diagonal(neighbours, False)
        ii, jj = np.where(neighbours)
        X, Y, Z = positions[ii].T
        U, V, W = (positions[jj] - positions[ii]).T
        C = rij[neighbours]
        C = C / C.min()

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d") if rank == 3 else fig.add_subplot()
    ax.set_aspect("equal")

    if rank == 3:
        if draw_bonds:
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

    else:
        ax.set_xlabel("x" + f" ({units})" if units is not None else "x")
        ax.set_ylabel("y" + f" ({units})" if units is not None else "y")

        if draw_bonds:
            ax.quiver(
                X,
                Y,
                U,
                V,
                linewidth=1,
                angles="xy",
                scale_units="xy",
                scale=2,
                headaxislength=0,
                headlength=0,
                color="gray",
                alpha=1 / C**3,
            )

        if colouring is not None:
            colours = ["C0", "C2", "C1", "C3", "C4", "C5", "C6"]  # green as 2nd
            for c, label in enumerate(set(colouring)):
                inds = [j == label for j in colouring]
                ax.scatter(x[inds], y[inds], c=colours[c], label=label, s=80)
            ax.legend()
        else:
            ax.scatter(x, y, c="g", s=80)

        if show_labels:
            for ind in range(qbits.nqbits):
                ax.text(x[ind], y[ind], s=qbits.labels[ind])
