"""
Quantum Simulation Environment.


This package is adapted from Atomic Simulation Environment (ASE).
"""

__all__ = [
    "calc",
    "cell",
    "draw",
    "lattices",
    "magnetic",
    "Operator",
    "Qbit",
    "Qbits",
    "Signal",
    "utils",
    "visualise",
]
import importlib.metadata

__version__ = importlib.metadata.version("qse")

from qse.qbit import Qbit
from qse.qbits import Qbits
from qse.signal import Signal
from qse.visualise import draw
from qse.gate_based import Operator

from qse import calc, lattices, magnetic, utils, visualise  # isort: skip
