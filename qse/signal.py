"""Definition of the Signal class."""

from typing import Union

import numpy as np

class Signal:
    """
    Signal class represents a 1D signal with values and duration.

    The class supports arithmetic operations with other signals and
    scalars. Specifically:

    - Adding a scalar to a Signal returns a new Signal with modified
    values and the same duration. Say W is a signal, then W + 3.0 returns
    a signal with same duration, and values W.values + 3.0.
    
    - Adding two Signal instances concatenates their values and sums
    their durations. So if w1, w2 are instantiation of Signal, then
    w = w1 + w2 gives signal with concatenated values, i.e.,
    w.values = [w1.values, w2.values], and added duration.
    w.duration = w1.duration + w2.duration

    Attributes
    ----------
    values : ndarray
        The values of the signal.
    duration : int
        Duration of the signal.

    Note
    ----
    Currently, the object gets created for multi-dim arrays as well.
    However, it should be used for 1D only, we haven't made it useful
    or consistent for multi-dim usage.
    """

    def __init__(self, values, duration=None) -> None:
        """Instantiates the Signal class.

        Parameters
        ----------
        values : array_like
            The signal values.
        duration : int, optional
            Duration of the signal, by default None
        """
        self.values = np.asarray(values)
        self._duration = len(self.values) if duration is None else int(duration)

    @property
    def duration(self) -> int:
        """
        Duration of the signal.

        Returns
        -------
        int
            The duration.
        """
        return self._duration

    @duration.setter
    def duration(self, value):
        self._duration = value

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
        return (
            (self.duration == other.duration)
        and (self.values == other.values).all()
        )

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
        if isinstance(other, Signal):
            return Signal(
                values=np.append(self.values, other.values),
                duration=self.duration + other.duration,
            )
        elif isinstance(other, (float, int)):
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
        if isinstance(other, Signal):
            self.values = np.append(self.values, other.values)
            self.duration += other.duration
        elif isinstance(other, (float, int)):
            self.values += other
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
            self.values *= other
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

    # TODO: Define interpolating scheme to resample points
    # if duration is changed externally.


