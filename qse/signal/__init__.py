"""

signal
======

This module contains `Signal` ans `Signals` classes
to define a a quantity as a function of time. It can
be used to express quantities like time dependent
external field coupled to Hamiltonian.


"""

__all__ = [
    "Signal",
    "Signals",
]

from .signal import Signal
from .signals import Signals
