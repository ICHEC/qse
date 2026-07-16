import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


def draw_scalar(qbits, scalar, radius=None, show_labels=False, units=None):
    """
    Visualize the positions of a set of qubits with some scalar quantity.

    Parameters
    ----------
    qbits: qse.Qbits
        The Qbits object.
    scalar: list | np.ndarray
        Must have the same length as the number of Qubits.
    radius: float | str
        A cutoff radius for visualizing bonds.
        Pass 'nearest' to set the radius to the smallest
        distance between the passed qubits.
        If no value is passed the bonds will not be visualized.
    show_labels: bool
        Whether to show the labels of the qubits.
        Defaults to False.
    units : str, optional
        The units of distance.
    """
    if len(scalar) != qbits.nqbits:
        raise Exception("The length of scalar must equal the number of Qubits.")

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
    ax = fig.add_subplot()

    ax.set_xlabel("x" + f" ({units})" if units is not None else "x")

    if qbits.dim == 2:
        x, y = qbits.positions.T
        ax.set_ylabel("y" + f" ({units})" if units is not None else "y")
        ax.set_aspect("equal")

    else:
        x = qbits.positions.T.flatten()
        y = np.zeros(qbits.nqbits)
        ax.set_yticks([])  # remove y-ticks for 1d plot.

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


    ax.scatter(x, y, c="k", s=110, alpha=0.5)
    scatter = ax.scatter(x, y, c=scalar, s=100, cmap=mpl.colormaps["Blues"])
    fig.colorbar(scatter)

    if show_labels:
        for ind in range(qbits.nqbits):
            ax.text(x[ind], y[ind], s=qbits.labels[ind])
