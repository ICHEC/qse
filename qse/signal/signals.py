"""Definition of the Signal class."""

import numpy as np

from .draw import draw
from .signal import Signal


class Signals:
    """
    The Signal class represents a collection of Signals.

    Parameters
    ----------
    signals : list[qse.Signal], optional
        The signals.

    Examples
    --------
    >>> x = qse.Signal([1], 10)
    >>> y = qse.Signal(np.linspace(0, 1, 5), 5)
    >>> qse.Signals([x, y])
    ... Total duration=15
    ...   Signal(duration=10, values=[1.])
    ...   Signal(duration=5, values=[0.   0.25 0.5  0.75 1.  ])

    We can also create Signals by addition, for example to create
    the same signal as above:

    >>> z = qse.Signals()
    >>> z += x
    >>> z += y
    """

    def __init__(self, signals=None):
        if signals is None:
            signals = []
        else:
            if not isinstance(signals, list) or not isinstance(signals[0], Signal):
                raise Exception("signals must be a list of Signal objects.")
        self._signals = signals

    def __add__(self, other):
        if isinstance(other, Signal):
            return Signals(self._signals + [other])
        else:
            raise TypeError(f"Unsupported operand type for +: {type(other)}")

    def __iadd__(self, other):
        if isinstance(other, Signal):
            self._signals += [other]
        else:
            raise TypeError(f"Unsupported operand type for +=: {type(other)}")
        return self

    def __getitem__(self, i):
        return self.signals[i]

    def __repr__(self) -> str:
        return f"Total duration={self.duration}\n" + "\n".join(
            [f"  Signal(duration={s.duration}, values={s.values})" for s in self]
        )

    def __len__(self):
        return len(self.signals)

    @property
    def signals(self):
        return self._signals

    @property
    def duration(self):
        return sum([signal.duration for signal in self])

    def to_pulser(self):
        """
        Convert to a Pulser Waveform.

        Returns
        -------
        pulser.waveforms.Waveform
            The waveform.
        """
        from pulser.waveforms import CompositeWaveform

        if len(self) == 1:
            return self[0].to_pulser()
        return CompositeWaveform(*[i.to_pulser() for i in self])

    def expand(self):
        """
        Get an array of length 'duration' whose entries are the signal values.

        Returns
        -------
        np.ndarray
            An array representing the signal.
        """
        return np.concatenate([i.expand() for i in self])

    def draw(self, time_units=None, signal_units=None):
        """
        Draw the signal.

        Parameters
        ----------
        time_units : str, optional
            The units of the duration.
        signal_units : str, optional
            The units of the signal.
        """
        return draw(self, time_units=time_units, signal_units=signal_units)
