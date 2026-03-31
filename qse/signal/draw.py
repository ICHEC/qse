import matplotlib.pyplot as plt


def draw(signal, time_units=None, signal_units=None):
    fig = plt.figure()
    y = signal.expand()
    plt.plot(y)
    plt.xlabel(f"Time ({time_units})" if time_units is not None else "Time")
    plt.ylabel(f"Signal ({signal_units})" if signal_units is not None else "Signal")
    return fig
