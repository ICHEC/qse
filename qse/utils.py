"""
Utility functions for a variety of usage within and outside of QSE.
"""

import numpy as np


def int2bin(x, width=32):
    """
    Converts an integer array to array of equivalent binary strings.

    Parameters
    ----------
    x : np.ndarray
        An array of integers.
    width : int
        The length of the binary strings.

    Returns
    -------
    np.ndarray
        The array of binary strings.

    Notes
    -----
    This function is equivalent to:

    >>> int2bin = np.vectorize(lambda x, width=16: np.binary_repr(x,width=width))

    However vectorize version is a bit slower compared to the one below.
    """
    out = np.fromiter(
        (np.binary_repr(i, width=width) for i in x), dtype=f"U{width}", count=x.shape[0]
    )
    return out
