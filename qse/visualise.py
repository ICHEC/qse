import matplotlib.pyplot as plt
import numpy as np


def draw(qbits, radius=None, draw_bonds=True, show_labels=False, colouring=None):
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
    show_labels: bool
        Whether to show the labels of the qubits.
        Defaults to False.
    colouring: str | list
        A set of integers used to assign different colors to each Qubit.
        This can be used to view different magnetic orderings.
        Must have the same length as the number of Qubits.
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
    print(f"cell_rank = {cell_rank}, position_rank = {position_rank}")
    positions = qbits.positions.copy()
    x, y, z = positions.T

    if draw_bonds:
        rij = qbits.get_all_distances()
        ib = np.logical_not(np.eye(qbits.nqbits, dtype=bool))
        rcut0 = rij[ib].min()

        if radius is None:
            rcut = rcut0
        else:
            if rcut0 > radius:
                rcut = radius
                draw_bonds = False
            else:
                rcut = radius

    if draw_bonds:
        nearest_neighbours = rij <= rcut
        np.fill_diagonal(nearest_neighbours, False)
        ii, jj = np.where(nearest_neighbours)
        X, Y, Z = positions[ii].T
        U, V, W = (positions[jj] - positions[ii]).T
        C = rij[nearest_neighbours]
        C = C / C.min()

    fig = plt.figure()

    if rank == 3:
        ax = fig.add_subplot(projection="3d") if rank == 3 else fig.add_subplot()
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

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
        ax = fig.add_subplot()
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

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
            for c, i in enumerate(set(colouring)):
                inds = [j == i for j in colouring]
                ax.scatter(x[inds], y[inds], c=f"C{c}", label=i)
            ax.legend()

        else:
            ax.scatter(x, y, c="g")
        if show_labels:
            for ind in range(qbits.nqbits):
                ax.text(x[ind], y[ind], s=qbits.labels[ind])
