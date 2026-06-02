import matplotlib.pyplot as plt
import numpy as np

from .colours import qse_palette

colors = qse_palette["colors"]
rads = qse_palette["rads"]


def draw_qbits(
    qbits,
    radius=None,
    show_labels=False,
    colouring=None,
    units=None,
    equal_aspect=True,
    alpha_min=0.0,
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
    alpha_min : float, optional
        Minimum alpha for bond opacity. Bond alphas are linearly rescaled
        from (alpha_min, 1), where 1 is the shortest bond and alpha_min
        is the longest. Defaults to 0.0.
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

    figsize = (10, 8) if qbits.dim == 3 else (6.4, 4.8)
    fig = plt.figure(figsize=figsize)
    projection = "3d" if qbits.dim == 3 else None
    ax = fig.add_subplot(projection=projection)
    if equal_aspect:
        ax.set_aspect("equal")

    if qbits.dim == 3:
        _draw_3d(
            qbits,
            draw_bonds,
            radius,
            rij,
            min_dist,
            alpha_min,
            units,
            colouring,
            show_labels,
            ax,
        )
    else:
        _draw_2d(
            qbits,
            draw_bonds,
            radius,
            rij,
            min_dist,
            alpha_min,
            units,
            colouring,
            show_labels,
            ax,
        )
    return fig


def _draw_3d(
    qbits,
    draw_bonds,
    radius,
    rij,
    min_dist,
    alpha_min,
    units,
    colouring,
    show_labels,
    ax,
):
    positions = qbits.positions

    ax.set_xlabel("x" + f" ({units})" if units is not None else "x")
    ax.set_ylabel("y" + f" ({units})" if units is not None else "y")
    ax.set_zlabel("z" + f" ({units})" if units is not None else "z")
    ax.figure.subplots_adjust(right=0.85)

    if draw_bonds:
        f_tol = 1.01  # fractional tolerance
        neighbours = rij <= radius * f_tol
        np.fill_diagonal(neighbours, False)
        ii, jj = np.where(neighbours)
        X, Y, Z = positions[ii].T
        U, V, W = (positions[jj] - positions[ii]).T
        alpha = alpha_min + (1 - alpha_min) * (min_dist / rij[neighbours]) ** 3

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

    if colouring is not None:
        inds0 = [j == 0 for j in colouring]
        inds1 = [j == 1 for j in colouring]
        for r, c in zip(rads, colors):
            ax.scatter(
                x[inds0],
                y[inds0],
                z[inds0],
                s=r**2,
                color=(0.1, c, 0.5),
                zorder=1,
                alpha=0.8,
            )
            ax.scatter(
                x[inds1],
                y[inds1],
                z[inds1],
                s=r**2,
                color=(c, 0.1, 0.5),
                zorder=1,
                alpha=0.8,
            )
    else:
        for r, c in zip(rads, colors):
            ax.scatter(x, y, z, s=r**2, color=(0.1, c, 0.5), zorder=1, alpha=0.8)

    if show_labels:
        for ind in range(qbits.nqbits):
            ax.text(x[ind], y[ind], z[ind], s=qbits.labels[ind])


def _draw_2d(
    qbits,
    draw_bonds,
    radius,
    rij,
    min_dist,
    alpha_min,
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
            alpha = alpha_min + (1 - alpha_min) * (min_dist / rij[i, j]) ** 3
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


def draw_3d_qbits_interactive(
    qbits, radius=None, show_labels=False, colouring=None, units=None, alpha_min=0.0
):
    """
    Visualize the positions of a set of qubits as an interactive 3D Plotly figure.

    This function is intended for 3D lattices only. For 2D or 1D qubits use
    :func:`draw_qbits` instead.

    Produces a rotatable, zoomable plot that renders as interactive HTML in
    Jupyter Book without requiring a live kernel.

    Parameters
    ----------
    qbits : qse.Qbits
        The Qbits object. Must have ``qbits.dim == 3``.
    radius : float | str, optional
        A cutoff radius for visualizing bonds.
        Pass 'nearest' to use the smallest distance between qubits.
        If None, bonds are not drawn.
    show_labels : bool, optional
        Whether to show qubit labels. Defaults to False.
    colouring : list, optional
        A list of 0s and 1s assigning each qubit to a sublattice.
        0 → green, 1 → red. Must have the same length as the number of qubits.
    units : str, optional
        The units of distance, shown on the axis labels.
    alpha_min : float, optional
        Minimum opacity for bonds. Defaults to 0.0.
    """
    import plotly.graph_objects as go

    if colouring is not None:
        if len(colouring) != qbits.nqbits:
            raise Exception("The length of colouring must equal the number of Qubits.")
        colouring = [int(i) for i in colouring]

    positions = qbits.positions
    x, y, z = positions.T

    def axis_label(name):
        return name + f" ({units})" if units is not None else name

    fig = go.Figure(
        layout=go.Layout(
            width=800,
            height=700,
            scene=dict(
                xaxis_title=axis_label("x"),
                yaxis_title=axis_label("y"),
                zaxis_title=axis_label("z"),
            ),
        )
    )

    # --- bonds ---
    draw_bonds = radius is not None
    rij = None
    if draw_bonds:
        rij = qbits.get_all_distances()
        min_dist = rij[np.logical_not(np.eye(qbits.nqbits, dtype=bool))].min()
        if radius == "nearest":
            radius = min_dist
        elif min_dist > radius:
            draw_bonds = False

    if draw_bonds:
        f_tol = 1.01
        x_lines, y_lines, z_lines = [], [], []
        for i in range(qbits.nqbits - 1):
            for j in range(i + 1, qbits.nqbits):
                if rij[i, j] <= radius * f_tol:
                    x_lines += [positions[i, 0], positions[j, 0], None]
                    y_lines += [positions[i, 1], positions[j, 1], None]
                    z_lines += [positions[i, 2], positions[j, 2], None]

        fig.add_trace(
            go.Scatter3d(
                x=x_lines,
                y=y_lines,
                z=z_lines,
                mode="lines",
                line=dict(color="gray", width=2),
                opacity=alpha_min + (1 - alpha_min) * 0.5,
                showlegend=False,
                hoverinfo="skip",
            )
        )

    # --- qubits ---
    mode = "markers+text" if show_labels else "markers"
    labels = list(qbits.labels) if show_labels else None

    def _rgb(r, g, b):
        return "rgb({},{},{})".format(int(r * 255), int(g * 255), int(b * 255))

    if colouring is not None:
        inds0 = np.array([c == 0 for c in colouring])
        inds1 = np.array([c == 1 for c in colouring])
        for rad, c in zip(rads, colors):
            fig.add_trace(
                go.Scatter3d(
                    x=x[inds0],
                    y=y[inds0],
                    z=z[inds0],
                    mode=mode,
                    marker=dict(size=rad / 2, color=_rgb(0.1, c, 0.5), opacity=0.8),
                    text=np.array(labels)[inds0] if labels else None,
                    showlegend=False,
                    hoverinfo="skip",
                )
            )
            fig.add_trace(
                go.Scatter3d(
                    x=x[inds1],
                    y=y[inds1],
                    z=z[inds1],
                    mode=mode,
                    marker=dict(size=rad / 2, color=_rgb(c, 0.1, 0.5), opacity=0.8),
                    text=np.array(labels)[inds1] if labels else None,
                    showlegend=False,
                    hoverinfo="skip",
                )
            )
    else:
        for rad, c in zip(rads, colors):
            fig.add_trace(
                go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode=mode,
                    marker=dict(size=rad / 2, color=_rgb(0.1, c, 0.5), opacity=0.8),
                    text=labels,
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    try:
        from IPython.display import HTML, display

        display(HTML(fig.to_html(full_html=False, include_plotlyjs="require")))
    except ImportError:
        fig.show()
