import numpy as np


class Cell:
    """
    A class to represent a lattice cell.

    Attributes
    ----------
    lattice_vectors : numpy.ndarray
        A 3x3 matrix representing the lattice vectors.
        Each row is a lattice vector.
        Defaults to all zeros.
    """

    def __init__(self, lattice_vectors=None):
        if lattice_vectors is None:
            lattice_vectors = np.zeros((3, 3))
        self.lattice_vectors = lattice_vectors

    def __repr__(self):
        return self.to_str()

    @property
    def lattice_vectors(self):
        return self._lattice_vectors

    @lattice_vectors.setter
    def lattice_vectors(self, value):
        value = np.array(value)
        if value.shape != (3, 3):
            raise Exception("The lattice vectors must be a 3x3 matrix.")
        self._lattice_vectors = value

    def rank(self):
        """
        Compute the rank of the cell.

        Returns
        -------
        int
            The rank of the cell.
        """
        return np.linalg.matrix_rank(self.lattice_vectors)

    def volume(self):
        """
        Compute the volume of the cell.

        Returns
        -------
        float
            The cell volume.
        """
        return np.linalg.det(self.lattice_vectors)

    def reciprocal(self):
        """
        Compute the reciprocal lattice vectors.

        Returns
        -------
        numpy.ndarray
            The reciprocal cell, calculated as 2π * inverse of the
            transposed cell matrix.
        """
        return 2 * np.pi * np.linalg.inv(self.lattice_vectors.T)

    def to_str(self):
        return "\n".join(
            [
                " ".join(["{:10.5f}".format(val) for val in row])
                for row in self.lattice_vectors
            ]
        )
