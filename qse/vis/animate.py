import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

viridis = mpl.colormaps["viridis"]


def animate_qbits(qbits, scalar, radius=None):
    """
    Animate the evolution of some scalar quantity.

    Parameters
    ----------
    qbits: qse.Qbits
        The Qbits object.
    scalar: list | np.ndarray
        The scalar quantity to define the animation.
        Must have the same length as the number of qubits.
    radius: float | str
        A cutoff radius for visualizing bonds.
        Pass 'nearest' to set the radius to the smallest
        distance between the passed qubits.
        If no value is passed the bonds will not be visualized.
    """
    if len(scalar) != qbits.nqbits:
        raise Exception("The length of scalar must equal the number of Qubits.")

    # Set up the figure and axis
    fig, ax = plt.subplots()

    if qbits.dim == 2:
        x, y = qbits.positions.T

    else:
        x = qbits.positions.T.flatten()
        y = np.zeros(qbits.nqbits)

    _draw_bonds(qbits, ax, radius, x, y)

    circles = []
    for q in qbits:
        circle = plt.Circle((q.x, q.y), 0.1, color=viridis(scalar[0, q.index]))
        ax.add_patch(circle)
        circles.append(circle)

    # Animation update function
    def update(frame):
        for q in qbits:
            circles[q.index].set_color(viridis(scalar[frame, q.index]))
        return circles

    ax.set_aspect("equal")
    ax.axis("off")  # Hide axes
    ax.set_xlim(x.min() - 0.5, x.max() + 0.5)
    ax.set_ylim(y.min() - 0.5, y.max() + 0.5)

    # Create the animation
    ani = FuncAnimation(
        fig,
        update,
        frames=len(scalar),
        interval=200,  # 50ms between frames
        blit=True,
    )
    return ani


def _draw_bonds(qbits, ax, radius, x, y):
    if radius is None:
        return

    rij = qbits.get_all_distances()
    min_dist = rij[np.logical_not(np.eye(qbits.nqbits, dtype=bool))].min()
    if radius == "nearest":
        radius = min_dist
    elif min_dist > radius:
        return

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
