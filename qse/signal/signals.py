"""Definition of the Signal class."""

from .signal import Signal


class Signals:
    """
    The Signal class represents a collection of Signals.

    Parameters
    ----------
    signals : list[qse.calc.Signal]
        The signals.
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
        from pulser.waveforms import CompositeWaveform

        if len(self) == 1:
            return self[0].to_pulser()
        return CompositeWaveform(*[i.to_pulser() for i in self])
