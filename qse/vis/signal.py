import matplotlib.pyplot as plt


def draw_signal(signal, time_units=None, signal_units=None, title=None):
    """
    Draw the signal.

    Parameters
    ----------
    signal : qse.Signal | qse.Signals | list[qse.Signal] | list[qse.Signals]
        The signal to be drawn. To view a single signal pass a qse.Signal or qse.Signals object.
        To view multiple signals pass a list of qse.Signal or qse.Signals objects.
    time_units : str, optional
        The units of the duration.
    signal_units : str, optional
        The units of the signal.
    title : str, optional
        A title for the plot.
    """
    fig = plt.figure()

    if not isinstance(signal, list):
        signal = [signal]

    for s in signal:
        plt.plot(s.expand(), label=s.name)
    
    plt.xlabel(f"Time ({time_units})" if time_units is not None else "Time")
    plt.ylabel(f"Signal ({signal_units})" if signal_units is not None else "Signal")

    if len([s.name for s in signal if s.name]) > 0:
        plt.legend()

    if title:
        plt.title(title)
    return fig
