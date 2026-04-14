import matplotlib.pyplot as plt


def draw_signal(signal, time_units=None, signal_units=None, title=None):
    """
    Draw the signal.

    Parameters
    ----------
    time_units : str, optional
        The units of the duration.
    signal_units : str, optional
        The units of the signal.
    title : str, optional
        A title for the plot.
    """
    fig = plt.figure()
    y = signal.expand()
    plt.plot(y)
    plt.xlabel(f"Time ({time_units})" if time_units is not None else "Time")
    plt.ylabel(f"Signal ({signal_units})" if signal_units is not None else "Signal")
    if title:
        plt.title(title)
    return fig
