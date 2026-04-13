r"""
QUBO
----

This module helps frame the QUBO (Quadratic Unconstrained
Binary Optimization) problem. The Qubo problem is finding
minima of the following function

.. math::

    F(x_1, x_2, ..., x_N) = \sum_{i, j} Q_{ij} x_i x_j


where :math:`x_i` are binary variables, :math:`i` runs over
the number of variables :math:`N` and the matrix :math:`Q_{ij}`
defines the QUBO problem.

Several interesting usecases map to QUBO, such as

- Ising model
- Set partitioning problem
- Portfolio optimization
- Traffic optimization

"""

# import matplotlib.pyplot as plt
import numpy as np

# from scipy.optimize import minimize
from scipy.spatial.distance import cdist  # , pdist, squareform

# import qse


class Qubo:
    """
    Qubo Class

    Parameters
    ----------
    N : int
        Problem size, for 1D lattice system,
        extent in one dimension. For 2D system
        this takes extent in first dimension
    N2 : int, optional
        If provides, the system is 2D, and N2
        takes the extent in second dimension, by default None
    C6 : float, optional
        Coefficient of interaction, by default 1.0
    Q : np.ndarray, optional
        The matrix for the QUBO problem, by default None
    """

    def __init__(self, N, N2=None, C6=1.0, Q=None) -> None:
        if N2 is None:
            self.N = N
            self.L1 = N
            self.L2 = 1
        else:
            self.L1 = N
            self.L2 = N2
            self.N = N * N2
        self.C6 = C6
        self.offdiagonal = ~np.eye(self.N, dtype=bool)

        if Q is None:
            self._Q = self.random_Q()
        else:
            self._Q = Q
        self._qrij = self.Q2rij()

    def nn1d_Q(self, diag=-1.0, offdiag=2.0):
        qq = np.eye(self.N, self.N)
        np.fill_diagonal(qq, diag)
        i, j = np.diag_indices(self.N)
        upper = (i, np.roll(j, -1))
        lower = (np.roll(j, -1), i)
        qq[upper] = offdiag
        qq[lower] = offdiag
        qq[0, -1] = 0
        qq[-1, 0] = 0
        return qq

    def radial1d_Q(self, diag=-1.0, offdiag=2.0, power=2) -> np.ndarray:
        qq = np.fromfunction(lambda i, j: abs(i - j) ** power, (self.N, self.N))
        np.fill_diagonal(qq, diag)
        qq = np.divide(1, qq, where=self.offdiagonal, out=None)
        return qq

    def radial2d_Q(self, diag=-1.0, offdiag=2.0, power=2) -> np.ndarray:
        qq = np.fromfunction(
            lambda i1, i2, j1, j2: abs(i1 - j1) ** power + abs(i2 - j2) ** power,
            (self.L1, self.L2, self.L1, self.L2),
        ).reshape(self.N, self.N)
        np.fill_diagonal(qq, diag)
        qq = np.divide(1, qq, where=self.offdiagonal, out=None)
        return qq

    def nn2d_Q(self, diag=-1.0, offdiag=2.0):
        qq = np.eye(self.N, self.N)
        np.fill_diagonal(qq, diag)
        indices = np.arange(self.L1 * self.L2).reshape(self.L1, self.L2)
        i1d = np.arange(self.L1)
        i2d = np.arange(self.L2)
        n1d = np.roll(i1d, -1)
        n2d = np.roll(i2d, -1)
        iis, nns = [], []
        for i1 in range(self.L1):
            for i2 in range(self.L2):
                ii = indices[i1, i2]
                j1, j2 = i1, n2d[i2]
                nns.append(indices[j1, j2])
                j1, j2 = n1d[i1], i2
                nns.append(indices[j1, j2])
                iis += [ii, ii]
        nns = np.array(nns)
        iis = np.array(iis)
        qq[(iis, nns)] = offdiag
        qq[(nns, iis)] = offdiag
        return qq

    def random_Q(self, diag=-1.0, offdiag=2.0) -> np.ndarray:
        qq = np.random.random((self.N, self.N))  # [ 0, 1]
        qq = qq + qq.T  # [ 0, 2]
        qq *= 0.5  # [ 0, 1]
        qq *= offdiag  # [-offdiag, offdiag]
        np.fill_diagonal(qq, diag)
        return qq

    def Q2rij(self) -> np.ndarray:
        rij = np.zeros((self.N, self.N), dtype=float)
        rij[:] = np.divide(
            self.C6 ** (1.0 / 6.0),
            np.where(self.offdiagonal, self.Q, 0.0),
            where=self.offdiagonal,
            out=None,
        )
        np.fill_diagonal(rij, 0)
        return rij

    def cost_function(self, x):
        """
        Cost function to minimize for mapping positions
        for a given Q

        Parameters
        ----------
        x : np.ndarray
            Flattened (original shape :math:`N \\times 3` as coordinates)
            array of length :math:`3N`
        """
        if self.N * 2 != x.shape[0] or len(x.shape) != 1:
            raise ValueError(f"x should of size {self.N * 2}, got {x.size}")
        Ri = np.reshape(x, (self.N, 2))
        rij = cdist(Ri, Ri)
        sum_square = ((self.qrij - rij) ** 2).sum()
        return sum_square

    @property
    def Q(self):
        return self._Q

    @Q.setter
    def Q(self, q):
        self._Q = q
        self._qrij = self.Q2rij()

    @property
    def qrij(self):
        return self._qrij
