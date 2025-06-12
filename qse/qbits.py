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

    def rotate(self, a, v, center=(0, 0, 0)):
        """
        Rotate qbits based on a vector and an angle, or two vectors.

        Parameters
        ----------
        a :
            Angle that the qbits is rotated around the vector 'v'. 'a'
            can also be a vector and then 'a' is rotated
            into 'v'.
        v :
            Vector to rotate the qbits around. Vectors can be given as
            strings: 'x', '-x', 'y', ... .
        center :
            The center is kept fixed under the rotation. Use 'COP' to
            fix the center of positions or 'COU' to fix the center of
            cell. Defaults to = (0, 0, 0).
        rotate_cell = False:
            If true the cell is also rotated.

        Examples
        --------
        Rotate 90 degrees around the z-axis, so that the x-axis is
        rotated into the y-axis:

        >>> qbits = Qbits()
        >>> qbits.rotate(90, 'z')
        >>> qbits.rotate(90, (0, 0, 1))
        >>> qbits.rotate(-90, '-z')
        >>> qbits.rotate('x', 'y')
        >>> qbits.rotate((1, 0, 0), (0, 1, 0))
        """

        if not isinstance(a, numbers.Real):
            a, v = v, a

        v = string2vector(v)

        normv = np.linalg.norm(v)

        if normv == 0.0:
            raise ZeroDivisionError("Cannot rotate: norm(v) == 0")

        if isinstance(a, numbers.Real):
            a *= pi / 180
            v /= normv
            c = cos(a)
            s = sin(a)
        else:
            v2 = string2vector(a)
            v /= normv
            normv2 = np.linalg.norm(v2)
            if normv2 == 0:
                raise ZeroDivisionError("Cannot rotate: norm(a) == 0")
            v2 /= np.linalg.norm(v2)
            c = np.dot(v, v2)
            v = np.cross(v, v2)
            s = np.linalg.norm(v)
            # In case *v* and *a* are parallel, np.cross(v, v2) vanish
            # and can't be used as a rotation axis. However, in this
            # case any rotation axis perpendicular to v2 will do.
            eps = 1e-7
            if s < eps:
                v = np.cross((0, 0, 1), v2)
                if np.linalg.norm(v) < eps:
                    v = np.cross((1, 0, 0), v2)
                assert np.linalg.norm(v) >= eps
            elif s > 0:
                v /= s

        center = self._centering_as_array(center)

        p = self.arrays["positions"] - center
        self.arrays["positions"][:] = (
            c * p - np.cross(p, s * v) + np.outer(np.dot(p, v), (1.0 - c) * v) + center
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

        Use mic=True to use the Minimum Image Convention and calculate
        the angle across periodic boundaries.
        """
        indices = np.array(indices)
        assert indices.shape[1] == 3

        a1s = self.positions[indices[:, 0]]
        a2s = self.positions[indices[:, 1]]
        a3s = self.positions[indices[:, 2]]

        v12 = a1s - a2s
        v32 = a3s - a2s

        return get_angles(v12, v32)

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
