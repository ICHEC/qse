"""Interface to different QSE calculators."""

__all__ = ["abc", "calculator", "signal"]
import qse.calc.abc
import qse.calc.calculator

from typing import Union

np = qse.np


class signal(object):
    """
    signal class to represent a a signal with a duration.
    It has two components: values, and duration.
    Instantiation of this class supports '+' in two ways.
    If w1, w2 are instantiation of signal, then
    - w1 + 3 gives a signal with values, w1.values + 3 and
    same duration.
    - w = w1 + w2 gives signal with concatenated values, i.e.,
    w.values = [w1.values, w2.values], and added duration.
    w.duration = w1.duration + w2.duration

    Note: Currently, the object gets created for multi-dim
    arrays as well. However, it should be used for 1D only,
    we haven't made it useful or consistent for multi-dim usage.
    """

    #
    # __slots__ = ('duration', 'values')
    def __init__(self, values, duration=None) -> None:
        self.values = np.asarray(values)
        self._duration = len(self.values) if duration is None else int(duration)

    #
    @property
    def duration(self):
        """time duration of signal"""
        return self._duration

    @duration.setter
    def duration(self, value):
        self._duration = value

    def __iter__(self):
        return iter(self.values)

    #
    def __getitem__(self, i):
        return self.values[i]

    #
    def __eq__(self, other) -> bool:
        return (self.duration == other.duration) and (self.values == other.values).all()

    #
    def __add__(self, other):
        if isinstance(other, signal):
            res = signal(
                values=np.append(self.values, other.values),
                duration=self.duration + other.duration,
            )
        else:
            if isinstance(other, Union[float, int]):
                res = signal(values=self.values + other, duration=self.duration)
            else:
                raise TypeError(f"wrong type for operand {type(other)}")
        return res

    #
    def __radd__(self, other):
        return self.__add__(other)

    #
    def __iadd__(self, other):
        if isinstance(other, signal):
            self.values = np.append(self.values, other.values)
            self.duration += other.duration
        else:
            if isinstance(other, Union[float, int]):
                self.values += other
            else:
                raise TypeError(f"wrong type for operand {type(other)}")
        return self

    #
    def __mul__(self, other):
        if isinstance(other, Union[float, int]):
            res = signal(values=self.values * other, duration=self.duration)
        else:
            raise TypeError(f"wrong type for operand {type(other)}")
        return res

    #
    def __rmul__(self, other):
        return self.__mul__(other)

    #
    def __imul__(self, other):
        if isinstance(other, Union[float, int]):
            self.values *= other
        else:
            raise TypeError(f"wrong type for operand {type(other)}")
        return self

    #
    def __repr__(self) -> str:
        return f"signal(duration={self.duration}, values={self.values})"

    # we need to define interpolating scheme to resample points if duration is changed externally.


#

from .pulser import Pulser
from .myqlm import Myqlm
