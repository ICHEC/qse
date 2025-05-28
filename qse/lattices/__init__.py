"""A collection of convenience functions for creating common lattices."""

__all__ = [
    "chain",
    "square",
    "triangular",
    "hexagonal",
    "kagome",
    "ring",
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
