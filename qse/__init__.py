"""
Quantum Simulation Environment.

This package is adapted from Atomic Simulation Environment (ASE).

Following are the major components of the module -

.. mermaid::

    mindmap
    root((QSE))
        qse.calc
        qse.lattices
        qse.magnetic
        qse.qbit
        qse.qbits
        qse.signal
        qse.utils


"""

__all__ = [
    "bar",
    "calc",
    "Cell",
    "draw",
    "lattices",
    "magnetic",
    "Qbit",
    "Qbits",
    "Operator",
    "Operators",
    "Signal",
    "Signals",
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
from qse.visualise import draw, bar

from qse import calc, lattices, magnetic, utils, visualise  # isort: skip
