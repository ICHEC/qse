import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


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
