"""Definition of the Signal class."""

import numpy as np


class Signal:
    """
    The Signal class represents a 1D signal with values and duration.

    Parameters
    ----------
    values : list | np.ndarray
        The values of the signal.
    duration : int
        Duration of the signal.
        Defaults to the length of the passed values.

    Examples
    --------
    One can create a signal by passing an array of values:

    >>> qse.Signal([1, 2, 3], 100)
    ... Signal(duration=100, values=[1. 2. 3.])

    Arithmetic operations with scalars is supported.
    Adding or multiplying a scalar to a Signal returns
    a new Signal with modified values and the same duration.
    For example:

    >>> signal = qse.Signal([1, 1])
    >>> signal * 3 + 0.5
    ... Signal(duration=2, values=[3.5 3.5])
    """

    def __init__(self, values, duration=None) -> None:
        self.values = np.asarray(values, dtype=float)
        self.duration = len(self.values) if duration is None else int(duration)

    def __iter__(self):
        """
        Iterate over the signal values.

        Returns
        -------
        iterator
            An iterator over the values.
        """
        return iter(self.values)

    def __getitem__(self, i):
        """
        Access individual signal value(s) by index or slice.

        Parameters
        ----------
        i : int or slice
            Index or slice to access.

        Returns
        -------
        scalar or ndarray
            The corresponding value(s).
        """
        return self.values[i]

    def __eq__(self, other) -> bool:
        """
        Check for equality with another signal.

        Parameters
        ----------
        other : Signal
            Signal to compare with.

        Returns
        -------
        bool
            True if equal in both values and duration.
        """
        return (self.duration == other.duration) and (self.values == other.values).all()

    def __add__(self, other):
        """
        Add another signal or scalar to this signal.

        Parameters
        ----------
        other : Signal or float or int
            Signal to concatenate or scalar to add elementwise.

        Returns
        -------
        Signal
            The resulting signal.

        Raises
        ------
        TypeError
            If the operand type is unsupported.
        """
        if isinstance(other, (float, int)):
            return Signal(values=self.values + other, duration=self.duration)
        else:
            raise TypeError(f"Unsupported operand type for +: {type(other)}")

    def __radd__(self, other):
        """
        Right addition operator.

        Returns
        -------
        Signal
            The resulting signal.
        """
        return self.__add__(other)

    def __iadd__(self, other):
        """
        In-place addition with another signal or scalar.

        Parameters
        ----------
        other : Signal or float or int
            Signal to concatenate or scalar to add elementwise.

        Returns
        -------
        Signal
            The updated signal.

        Raises
        ------
        TypeError
            If the operand type is unsupported.
        """
        if isinstance(other, (float, int)):
            self.values = self.values + other
        else:
            raise TypeError(f"Unsupported operand type for +=: {type(other)}")
        return self

    def __mul__(self, other):
        """
        Multiply signal values by a scalar.

        Parameters
        ----------
        other : float or int
            Scalar multiplier.

        Returns
        -------
        Signal
            The resulting signal.

        Raises
        ------
        TypeError
            If the operand type is unsupported.
        """
        if isinstance(other, (float, int)):
            return Signal(values=self.values * other, duration=self.duration)
        else:
            raise TypeError(f"Unsupported operand type for *: {type(other)}")

    def __rmul__(self, other):
        """
        Right multiplication by a scalar.

        Returns
        -------
        Signal
            The resulting signal.
        """
        return self.__mul__(other)

    def __imul__(self, other):
        """
        In-place multiplication by a scalar.

        Parameters
        ----------
        other : float or int
            Scalar multiplier.

        Returns
        -------
        Signal
            The updated signal.

        Raises
        ------
        TypeError
            If the operand type is unsupported.
        """
        if isinstance(other, (float, int)):
            self.values = self.values * other
        else:
            raise TypeError(f"Unsupported operand type for *=: {type(other)}")
        return self

    def __repr__(self) -> str:
        """
        Return a string representation of the signal.

        Returns
        -------
        str
            String representation.
        """
        return f"Signal(duration={self.duration}, values={self.values})"

    def __len__(self):
        return len(self.values)

    @property
    def duration(self):
        return self._duration

    @duration.setter
    def duration(self, new_duration):
        if not isinstance(new_duration, int):
            raise ValueError("The duration must be an ints")

        if new_duration % len(self) != 0:
            raise ValueError("The number of values must divide the duration.")

        self._duration = new_duration

    def time_per_value(self):
        return self.duration // len(self)

    def concatenate(self):
        return np.concatenate([[i] * self.time_per_value() for i in self.values])

    def to_pulser(self):
        from pulser.waveforms import ConstantWaveform, CustomWaveform

        if len(self) == 1:
            return ConstantWaveform(duration=self.duration, value=self.values[0])
        return CustomWaveform(self.concatenate())


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
