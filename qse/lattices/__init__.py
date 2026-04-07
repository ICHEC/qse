"""
A collection of convenience functions for creating common lattices.
Qbits class is very flexible to generate arbitrary structures. For
convenience, helper functions to generate several common types of
lattices is defined here. See the examples below -


.. jupyter-execute::

    import qse
    qsqr = qse.lattices.square(
        lattice_spacing=1.0,
        repeats_x=4, repeats_y=4)
    qsqr.draw(radius=2)


"""

__all__ = [
    "chain",
    "hexagonal",
    "kagome",
    "ring",
    "square",
    "torus",
    "triangular",
]

from .lattices import (
    chain,
    hexagonal,
    kagome,
    ring,
    square,
    torus,
    triangular,
)
