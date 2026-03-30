"""
Quantum Simulation Environment.


This package is adapted from Atomic Simulation Environment (ASE).

"""

__all__ = [
    "calc",
    "Cell",
    "draw",
    "lattices",
    "magnetic",
    "Qbit",
    "Qbits",
    "Operator",
    "Operators",
    "utils",
    "visualise",
]
import importlib.metadata

__version__ = importlib.metadata.version("qse")

from qse.cell import Cell
from qse.operator import Operator, Operators
from qse.qbit import Qbit
from qse.qbits import Qbits
from qse.signal import Signal, Signals
from qse.visualise import draw

from qse import calc, lattices, magnetic, utils, visualise  # isort: skip
