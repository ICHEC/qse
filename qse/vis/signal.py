import matplotlib.pyplot as plt

from .colours import qse_green, qse_red


def draw_signal(signal, time_units=None, signal_units=None, title=None):
    """
    Draw the signal.

    Parameters
    ----------
    signal : qse.Signal | qse.Signals
        The signal to be drawn.
    time_units : str, optional
        The units of the duration.
    signal_units : str, optional
        The units of the signal.
    title : str, optional
        A title for the plot.
    """
    fig = plt.figure()
    plt.plot(signal.expand(), c=qse_green)
    plt.xlabel(f"Time ({time_units})" if time_units is not None else "Time")
    plt.ylabel(f"Signal ({signal_units})" if signal_units is not None else "Signal")
    if title:
        plt.title(title)
    return fig


def draw_amp_and_det(
    amplitude, detuning, time_units=None, signal_units=None, title=None
):
    """
    Draw a signal amplitude together with a signal detuning.

    Parameters
    ----------
    amplitude : qse.Signal | qse.Signals
        The amplitude to be drawn.
    detuning : qse.Signal | qse.Signals
        The detuning to be drawn.
    time_units : str, optional
        The units of the duration.
    signal_units : str, optional
        The units of the signal.
    title : str, optional
        A title for the plot.
    """
    fig = plt.figure()

    plt.plot(amplitude.expand(), c=qse_green, label="Amplitude")
    plt.plot(detuning.expand(), c=qse_red, label="Detuning")

    plt.xlabel(f"Time ({time_units})" if time_units is not None else "Time")
    plt.ylabel(f"Signal ({signal_units})" if signal_units is not None else "Signal")

    if title:
        plt.title(title)

    plt.axhline(y=0, ls="--", c="k")
    plt.legend()

    return fig
