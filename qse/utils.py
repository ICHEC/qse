"""
Utility functions for a variety of usage within and outside of QSE.

More to elaborate later on.
"""

import numpy as np

import qse


def int2bin(x, width=32):
    """
    converts an integer array to array of equivalent binary strings.
    Equivalent to:
        int2bin = np.vectorize(lambda x, width=16: np.binary_repr(x,width=width))
    However vectorize version is a bit slower compared to the one below.
    """
    out = np.fromiter(
        (np.binary_repr(i, width=width) for i in x), dtype=f"U{width}", count=x.shape[0]
    )
    return out


a0 = 1.0
zerovec = np.zeros((1, 3))
sqr_cell = qse.cell.Cell([[a0, 0, 0], [0, a0, 0], [0, 0, 0]])
tri_cell = qse.cell.Cell([[a0, 0, 0], [a0 / 2, np.sqrt(3) * a0 / 2, 0], [0, 0, 0]])
hex_cell = qse.cell.Cell(
    [
        [3 * a0 / 2, np.sqrt(3) * a0 / 2, 0],
        [3 * a0 / 2, -np.sqrt(3) * a0 / 2, 0],
        [0, 0, 0],
    ]
)
kag_cell = qse.cell.Cell([[2 * a0, 0, 0], [a0, np.sqrt(3) * a0, 0], [0, 0, 0]])


def linear(a=a0, N=6):
    """
    Generate qse.Qbits object in linear chain geometry.
    """
    cell = qse.cell.Cell([[a0, 0, 0], [0, 0, 0], [0, 0, 0]])
    cell = (a / a0) * cell
    unit = qse.Qbits(positions=zerovec)
    unit.cell = cell
    linlattice = unit.repeat((N, 1, 1))
    return linlattice


def squarelattice(a=a0, N1=2, N2=2):
    """
    Generate qse.Qbits object in square lattice geometry.
    """
    sqcell = (a / a0) * sqr_cell
    unit = qse.Qbits(positions=zerovec)
    unit.cell = sqcell
    square_lattice = unit.repeat((N1, N2, 1))
    return square_lattice


def triangularlattice(a=a0, N1=2, N2=2):
    """
    Generate qse.Qbits object in triangular lattice geometry.
    """
    tricell = (a / a0) * tri_cell
    unit = qse.Qbits(positions=zerovec)
    unit.cell = tricell
    triangle_lattice = unit.repeat((N1, N2, 1))
    return triangle_lattice


def hexagonlattice(a=a0, N1=2, N2=2):
    """
    Generate qse.Qbits object in hexagonal lattice geometry.
    """
    hexcell = (a / a0) * hex_cell
    unit = qse.Qbits(positions=np.array([[0, 0, 0], [a, 0, 0]]))
    unit.cell = hexcell
    hexagon_lattice = unit.repeat((N1, N2, 1))
    return hexagon_lattice


def kagomelattice(a=a0, N1=2, N2=2):
    """
    Generate qse.Qbits object in kagome lattice geometry.
    """
    kagcell = (a / a0) * kag_cell
    unit = qse.Qbits(
        positions=np.array([[0, 0, 0], [a, 0, 0], [a / 2, np.sqrt(3) * a / 2, 0]])
    )
    unit.cell = kagcell
    kagome_lattice = unit.repeat((N1, N2, 1))
    return kagome_lattice
