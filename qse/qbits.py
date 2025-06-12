import numpy as np


class Qbits:
    def __init__(self, positions):
        self._positions = positions

    @property
    def positions(self):
        return self._positions.copy()

    @positions.setter
    def positions(self, positions):
        self._positions = positions

    def __len__(self):
        return len(self.positions)

    def __add__(self, other):
        qbits = self.copy()
        qbits += other
        return qbits

    def __iadd__(self, other):
        self.extend(other)
        return self

    def copy(self):
        """Return a copy."""
        return self.__class__(positions=self.positions)

    def extend(self, other):
        """Extend qbits object by appending qbits from *other*."""
        assert isinstance(other, Qbits)

        self.positions = np.concatenate([self.positions, other.positions])

    def append(self, qbit):
        """Append qbit to end."""
        self.extend(self.__class__([qbit]))

    def __delitem__(self, indices):
        """
        Delete a subset of the qbits.

        Parameters
        ----------
        indices : int | list
            The indices to be deleted.
        """
        if isinstance(indices, list) and len(indices) > 0:
            # Make sure a list of booleans will work correctly and not be
            # interpreted at 0 and 1 indices.
            indices = np.array(indices)

        mask = np.ones(len(self), bool)
        mask[indices] = False

        self.positions = self.positions[mask]

    def translate(self, displacement):
        """
        Translate qbit positions.

        Parameters
        ----------
        displacement : float | np.ndarray
            The displacement argument can be a float an xyz vector or an
            nx3 array (where n is the number of qbits).
        """
        self.positions += np.array(displacement)

    def get_centroid(self):
        """
        Get the centroid of the positions.

        Notes
        -----
        For a set of $k$ positions $\textbf{x}_1, \textbf{x}_2, ..., \textbf{x}_k$
        the centroid is given by
        $\frac{\textbf{x}_1 + \textbf{x}_2 + ... + \textbf{x}_k}{k}.$
        """
        return self.positions.mean(0)

    def set_centroid(self, centroid):
        """
        Set the centroid of the positions.

        Parameters
        ----------
        centroid : float | np.ndarray
            The new centroid. Can be a float or a xyz vector

        Notes
        -----
        For a set of $k$ positions $\textbf{x}_1, \textbf{x}_2, ..., \textbf{x}_k$
        the centroid is given by
        $\frac{\textbf{x}_1 + \textbf{x}_2 + ... + \textbf{x}_k}{k}.$
        """
        difference = centroid - self.get_centroid()
        self.positions += difference

    def rotate(self, angle: float, vector="z", center=(0, 0, 0)):
        """
        Rotate qbits based on an angle and vector.

        Parameters
        ----------
        angle : float
            Angle that the qbits is rotated around the vector.
        vector :
            Vector to rotate the qbits around. Vectors can be given as
            strings: 'x', '-x', 'y', ... . Defaults to 'z'.
        center :
            The center is kept fixed under the rotation. Use 'COP' to
            fix the center of positions or 'COU' to fix the center of
            cell. Defaults to = (0, 0, 0).

        Examples
        --------
        Rotate 90 degrees around the z-axis, so that the x-axis is
        rotated into the y-axis:

        >>> qbits = Qbits()
        >>> qbits.rotate(90, 'z')
        >>> qbits.rotate(90, (0, 0, 1))
        >>> qbits.rotate(-90, '-z')
        """
        if isinstance(vector, str):
            vector_dict = {
                "x": [1.0, 0.0, 0.0],
                "y": [0.0, 1.0, 0.0],
                "z": [0.0, 0.0, 1.0],
            }
            vector = vector_dict[vector]

        vector = np.array(vector)
        vector = _create_unit_vector(vector)

        angle *= np.pi / 180
        c = np.cos(angle)
        s = np.sin(angle)

        center = self._centering_as_array(center)

        p = self.positions - center
        self.positions = (
            c * p
            - np.cross(p, s * vector)
            + np.outer(np.dot(p, vector), (1.0 - c) * vector)
            + center
        )

    def _centering_as_array(self, center):
        if isinstance(center, str):
            if center.lower() == "centre":
                return self.get_centroid()
            else:
                raise ValueError("Cannot interpret center")
        return np.array(center, float)

    def get_angle(self, index_1: int, index_2: int, index_3: int):
        """
        Get the angle in degress formed by three qbits.

        Parameters
        ----------
        index_1 : int
            The index of the first qubit.
        index_2 : int
            The index of the second qubit.
        index_3 : int
            The index of the third qubit.

        Notes
        -----
        Let x1, x2, x3 be the vectors describing the positions of the three
        qubits. Then we calcule the angle between x1-x2 and x3-x2.
        """
        return self.get_angles([[index_1, index_2, index_3]])[0]

    def get_angles(self, indices):
        """
        Get angle formed by three qbits for multiple groupings.

        Calculate angle in degrees between vectors between qbits a2->a1
        and a2->a3, where a1, a2, and a3 are in each row of indices.
        """
        indices = np.array(indices)
        assert indices.shape[1] == 3

        a1s = self.positions[indices[:, 0]]
        a2s = self.positions[indices[:, 1]]
        a3s = self.positions[indices[:, 2]]

        v12 = _create_unit_vector(a1s - a2s)
        v32 = _create_unit_vector(a3s - a2s)

        return np.arccos(np.dot(v32, v12))

    def get_distance(self, a0, a1, vector=False):
        """
        Return distance between two qbits.

        Use mic=True to use the Minimum Image Convention.
        vector=True gives the distance vector (from a0 to a1).
        """
        return self.get_distances(a0, [a1], vector=vector)[0]

    def get_distances(self, a, indices, vector=False):
        """
        Return distances of qbit No.i with a list of qbits.

        vector=True gives the distance vector (from a to self[indices]).
        """
        p1 = self.positions[a]
        p2 = self.positions[indices]

        if vector:
            return np.array([p1 - p for p in p2])
        return np.array([np.linalg.norm(p1 - p) for p in p2])

    def get_all_distances(self):
        """
        Return distances of all of the qbits with all of the qbits.
        """
        distances = np.zeros((len(self),) * 2)
        for i in range(len(self) - 1):
            for j in range(i + 1, len(self)):
                distances[i, j] = np.linalg.norm(self.positions[i] - self.positions[j])
                distances[j, i] = distances[i, j]
        return distances


def _create_unit_vector(vector):
    normv = np.linalg.norm(vector)

    if normv == 0.0:
        raise ZeroDivisionError("Vector must have nonzero length.")

    return vector / normv
