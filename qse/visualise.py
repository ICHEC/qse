import matplotlib.pyplot as plt
import numpy as np


def draw(qbits, ax=None, radius=None):
    """
    ...

    Parameters
    ----------
    qbits: Qbits
        ...
    ax: ...
        ...
    radius: float
        The radius

    Returns
    -------
    plt
    """
    crank = qbits.cell.rank
    prank = np.linalg.matrix_rank(qbits.positions)
    # if crank is more than prank, it means that a higher dimensional cell is present, and actual
    # structure requires repettition of cells to clearly visualize. if prank is more than crank,
    # it means the repettion happens in lower dimension of a local geometry which visualized in
    # higher dimension. Either way, we need to see things in higher dimension.
    # TODO: Current version has a bug, in which some of the bonds aren't drawn in rectangular geometry.
    # Need to be looked at.

    rank = max(crank, prank)
    print(f"crank = {crank}, prank = {prank}")
    positions = qbits.positions.copy()
    x, y, z = positions.T

    rij = qbits.get_all_distances()
    ib = np.logical_not(np.eye(qbits.nqbits, dtype=bool))
    rcut0 = rij[ib].min()

    if radius is None:
        draw_bond = True
        rcut = rcut0
    else:
        if rcut0 > radius:
            rcut = radius
            draw_bond = False
        else:
            rcut = radius
            draw_bond = True

    print(f"rcut is {rcut}")
    nns = rij <= rcut
    np.fill_diagonal(nns, False)
    ii, jj = np.where(nns)
    if draw_bond:
        X, Y, Z = positions[ii].T
        U, V, W = (positions[jj] - positions[ii]).T
        C = rij[nns]
        C = C / C.min()
    #
    # print(x.shape, y.shape, z.shape, X.shape, Y.shape, Z.shape, U.shape, V.shape, W.shape)
    fig = plt.figure() if rank == 3 or rank == 2 else plt.figure(figsize=(12, 3))

    ax = fig.add_subplot(projection="3d") if rank == 3 else fig.add_subplot()
    # set the aspect ratio to 1, so that the ratio is actually
    # proportional to data size in x, y, z direction.
    if rank == 1:
        ax.set_aspect(len(qbits))
    else:
        ax.set_aspect("equal")
    # ax.set_xlabel('x', labelpad=-12)
    # ax.set_ylabel('y', labelpad=-12)
    ax.set_xticks([])
    ax.set_yticks([])
    if rank == 3:  # ax.set_zlabel('z', labelpad=-12)
        ax.set_zticks([])

    if rank == 3:
        if draw_bond:
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
    if rank == 2:
        if draw_bond:
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
        ax.plot(x, y, "o", color="green")
    if rank == 1:
        # y = np.zeros(x.shape)
        if draw_bond:
            ax.quiver(
                X,
                Y,
                U,
                V,
                linewidth=1,
                angles="xy",
                scale_units="xy",
                scale=1,
                headaxislength=0,
                headlength=0,
                color="gray",
                alpha=1 / C**3,
            )
        ax.plot(x, y, "o", ms=9, color="red")
