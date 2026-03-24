"""
This module defines the central object in the QSE package: the Qbits object.
"""

import numbers
import warnings
from math import cos, sin

import numpy as np

from qse.cell import Cell
from qse.operator import Operator, Operators
from qse.qbit import Qbit
from qse.visualise import draw as _draw


class Qbits:
    """
    The Qbits object can represent an isolated molecule, or a
    periodically repeated structure.  It has a unit cell and
    there may be periodic boundary conditions along any of the three
    unit cell axes. Information about the qbits (qubit state and
    position) is stored in ndarrays.

    Parameters
    ----------
    positions: list of xyz-positions
        Qubit positions.  Anything that can be converted to an
        ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
        ...].
    labels: list of str
        A list of strings corresponding to a label for each qubit.
    states: list of 2-length arrays.
        State of each qubit.
    cell: qse.Cell | np.ndarray
        Unit cell vectors.
        Either pass a matrix where each row corresponds to a lattice vector
        or a qse.Cell object.
        Default value is all zeros.
    pbc: one or three bool
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        value: False.
    calculator: calculator object
        Used to attach a calculator for doing computation.

    Examples
    --------
    Empty Qbits object:

    >>> qs = qse.Qbits()

    These are equivalent:

    >>> a = qse.Qbits(
    ...     positions=np.array([(0, 0, 0), (0, 0, 2)])
    ...     labels=['qb1', 'qb2'],
    ... )
    >>> a = qse.Qbits.from_qbit_list(
    ...     [Qbit('qb1', position=(0, 0, 0)), Qbit('qb2', position=(0, 0, 2))]
    ... )

    >>> xd = np.array(
    ...    [[0, 0, 0],
    ...     [0.5, 0.5, 0.5]]
    ... )
    >>> qdim = qse.Qbits(xd)
    >>> qdim.cell = np.eye(3)
    >>> qdim.pbc = True
    >>> qlat = qdim.repeat([3,3,3])

    The qdim will have shape = (1,1,1) and qlat will have shape = (3, 3, 3)

    Notes
    -----
    In order to do computation, a calculator object has to attached
    to the qbits object.
    """

    def __init__(
        self,
        positions=None,
        labels=None,
        states=None,
        cell=None,
        pbc=None,
        calculator=None,
    ):
        if labels is None:
            nqbits = 0 if positions is None else len(positions)
        else:
            if not isinstance(labels, list):
                raise Exception("'labels' must be a list.")
            nqbits = len(labels)

            if (positions is not None) and (len(positions) != nqbits):
                raise Exception(
                    "Both 'positions' and 'labels' must have the same length."
                )

        self.arrays = {}

        # labels
        if labels is None:
            labels = [str(i) for i in range(nqbits)]
        # We allow for labels up to length 12.
        self.new_array("labels", labels, "<U12")

        # positions
        if positions is None:
            positions = np.zeros((nqbits, 3))

        self.new_array("positions", positions, float)

        # cell
        if cell is None:
            self._cell = None
        else:
            self.cell = cell

        # states
        if states is None:
            states = np.zeros((positions.shape[0], 2), dtype=complex)
            states[:, 0] = 1
        self.new_array("states", states, complex, (2,))

        # shape
        self.shape = (1,) * self.dim

        # pbc
        self._pbc = np.zeros(3, bool)
        if pbc is None:
            pbc = False
        self.set_pbc(pbc)

        # calculator
        self.calc = calculator

    @classmethod
    def from_qbit_list(self, qbit_list):
        # Get data from a list or tuple of Qbit objects:
        data = {
            f"{name}s": [qbit.get_raw(name) for qbit in qbit_list]
            for name in ["label", "state", "position"]
        }
        return Qbits(**data)

    @property
    def calc(self):
        """Calculator object."""
        return self._calc

    @calc.setter
    def calc(self, calc):
        self._calc = calc
        if hasattr(calc, "set_qbits"):
            calc.set_qbits(self)

    @property
    def positions(self):
        return self.arrays["positions"]

    @positions.setter
    def positions(self, pos):
        self.arrays["positions"][:] = pos

    @property
    def dim(self):
        return self.positions.shape[1]

    @property
    def nqbits(self):
        return len(self)

    @property
    def cell(self):
        return self._cell

    @cell.setter
    def cell(self, new_cell):
        if not isinstance(new_cell, Cell):
            new_cell = Cell(new_cell)

        if new_cell.dim != self.dim:
            raise Exception(
                "The dimension of the cell must match the dimension of the qubits."
            )

        self._cell = new_cell

    @property
    def pbc(self):
        """Reference to pbc-flags for in-place manipulations."""
        return self._pbc

    @pbc.setter
    def pbc(self, pbc):
        self._pbc[:] = pbc

    def set_pbc(self, pbc):
        """Set periodic boundary condition flags."""
        self.pbc = pbc

    def get_pbc(self):
        """Get periodic boundary condition flags."""
        return self.pbc.copy()

    def new_array(self, name, a, dtype=None, shape=None):
        """Add new array.

        If *shape* is not *None*, the shape of *a* will be checked."""

        if dtype is not None:
            a = np.array(a, dtype, order="C")
            if len(a) == 0 and shape is not None:
                a.shape = (-1,) + shape
        else:
            if not a.flags["C_CONTIGUOUS"]:
                a = np.ascontiguousarray(a)
            else:
                a = a.copy()

        if name in self.arrays:
            raise RuntimeError("Array {} already present".format(name))

        for b in self.arrays.values():
            if len(a) != len(b):
                raise ValueError(
                    'Array "%s" has wrong length: %d != %d.' % (name, len(a), len(b))
                )
            break

        # We allow positions to be 1D, 2D or 3D
        if name == "positions":
            if len(a.shape) != 2 or a.shape[1] not in [1, 2, 3]:
                raise ValueError(
                    "Positions must have shape (n_q, d) "
                    "where n_q is the number of qbits and the "
                    "dimnension d is 1, 2 or 3."
                )

        if shape is not None and a.shape[1:] != shape:
            raise ValueError(
                'Array "%s" has wrong shape %s != %s.'
                % (name, a.shape, (a.shape[0:1] + shape))
            )

        self.arrays[name] = a

    def get_array(self, name, copy=True):
        """Get an array.

        Returns a copy unless the optional argument copy is false.
        """
        if copy:
            return self.arrays[name].copy()
        else:
            return self.arrays[name]

    def set_array(self, name, a, dtype=None, shape=None):
        """
        Update array.

        If *shape* is not *None*, the shape of *a* will be checked.
        If *a* is *None*, then the array is deleted."""

        b = self.arrays.get(name)
        if b is None:
            if a is not None:
                self.new_array(name, a, dtype, shape)
        else:
            if a is None:
                del self.arrays[name]
            else:
                a = np.asarray(a)
                if a.shape != b.shape:
                    raise ValueError(
                        'Array "%s" has wrong shape %s != %s.'
                        % (name, a.shape, b.shape)
                    )
                b[:] = a

    def copy(self):
        """Return a copy."""
        qbits = self.__class__(cell=self.cell, pbc=self.pbc)

        qbits.arrays = {}
        for name, a in self.arrays.items():
            qbits.arrays[name] = a.copy()

        qbits.shape = self.shape  # this was necessary, and took long time to realise!
        return qbits

    def todict(self):
        """For basic JSON (non-database) support."""
        d = {}
        d["labels"] = self.arrays["labels"]
        d["positions"] = self.arrays["positions"]
        d["states"] = self.arrays["states"]
        d["cell"] = self.cell
        d["pbc"] = self.pbc

        return d

    @classmethod
    def fromdict(cls, dct):
        """Rebuild qbits object from dictionary representation (todict)."""
        dct = dct.copy()
        kw = {}
        for name in ["labels", "positions", "states", "cell", "pbc"]:
            kw[name] = dct.pop(name)

        qbits = cls(**kw)
        nqbits = len(qbits)

        # Some arrays are named differently from the qbits __init__ keywords.
        # Also, there may be custom arrays.  Hence we set them directly:
        for name, arr in dct.items():
            assert len(arr) == nqbits, name
            assert isinstance(arr, np.ndarray)
            qbits.arrays[name] = arr
        return qbits

    def __len__(self):
        return len(self.arrays["positions"])

    def __repr__(self):
        """We use pymatgen type of style to print Qbits class"""
        tokens = []
        insep = " "
        N = len(self)
        if N < 33:
            for i in self:
                tokens.append("{0}".format(i) + ",\n")
        else:
            for i in self[:3]:
                tokens.append("{0}".format(i) + ",\n")
            tokens.append("...\n")
            for i in self[-3:]:
                tokens.append("{0}".format(i) + ",\n")

        if self.pbc.any() and not self.pbc.all():
            txt = "pbc={0}".format(self.pbc.tolist())
        else:
            txt = "pbc={0}".format(self.pbc[0])
        tokens.append(txt + ",\n")

        if self.cell is not None:
            tokens.append("cell=\n{0}".format(self.cell.to_str()) + ",\n")

        if self._calc is not None:
            tokens.append("calculator={0}".format(self._calc.__class__.__name__))
        txt = "{0}(\n{1}{2})".format(
            self.__class__.__name__, insep, insep.join(tokens) + "..."
        )
        return txt

    def __add__(self, other):
        qbits = self.copy()
        qbits += other
        return qbits

    def extend(self, other):
        """Extend qbits object by appending qbits from *other*."""
        if isinstance(other, Qbit):
            other = self.from_qbit_list([other])

        if self.dim != other.dim:
            raise Exception("Cannot add systems of differing dimensions.")

        n1 = len(self)
        n2 = len(other)

        for name, a1 in self.arrays.items():
            a = np.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = np.empty((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            a[:n1] = 0

            self.set_array(name, a)

    def __iadd__(self, other):
        self.extend(other)
        return self

    def append(self, qbit):
        """Append qbit to end."""
        self.extend(self.__class__([qbit]))

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __getitem__(self, indices):
        """
        Return a subset of the qbits.

        Parameters
        ----------
        indices : int | list | slice
            The indices to be returned.

        Returns
        -------
        Qbit | Qbits.
            If indices is a scalar a Qbit object is returned. If indices
            is a list or a slice, a Qbits object with the same cell, pbc, and
            other associated info as the original Qbits object is returned.
        """

        if isinstance(indices, numbers.Integral):
            nqbits = len(self)
            if indices < -nqbits or indices >= nqbits:
                raise IndexError("Index out of range.")
            return Qbit(qbits=self, index=indices)

        if not isinstance(indices, slice):
            indices = np.array(indices)
            # if indices is a mask.
            if indices.dtype == bool:
                if len(indices) != len(self):
                    raise IndexError(
                        f"Length of mask {len(indices)} must equal "
                        f"number of qbits {len(self)}"
                    )
                indices = np.arange(len(self))[indices]

        qbits = self.__class__(cell=self.cell, pbc=self.pbc)

        qbits.arrays = {}
        for name, a in self.arrays.items():
            qbits.arrays[name] = a[indices].copy()

        return qbits

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
        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]

    def pop(self, i=-1):
        """Remove and return qbit at index *i* (default last)."""
        qbit = self[i]
        qbit.cut_reference_to_qbits()
        del self[i]
        return qbit

    def __imul__(self, m):
        """In-place repeat of qbits."""
        if isinstance(m, int):
            m = (m,) * self.dim

        for x, vec in zip(m, self.cell.lattice_vectors):
            if x != 1 and not vec.any():
                raise ValueError("Cannot repeat along undefined lattice " "vector")

        fac = np.prod(m)
        n = len(self)

        for name, a in self.arrays.items():
            self.arrays[name] = np.tile(a, (fac,) + (1,) * (len(a.shape) - 1))

        positions = self.arrays["positions"]

        i0 = 0

        if self.dim == 1:
            for m0 in range(m[0]):
                i1 = i0 + n
                positions[i0:i1] += np.dot((m0), self.cell.lattice_vectors)
                i0 = i1

        elif self.dim == 2:
            for m0 in range(m[0]):
                for m1 in range(m[1]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1), self.cell.lattice_vectors)
                    i0 = i1

        elif self.dim == 3:
            for m0 in range(m[0]):
                for m1 in range(m[1]):
                    for m2 in range(m[2]):
                        i1 = i0 + n
                        positions[i0:i1] += np.dot(
                            (m0, m1, m2), self.cell.lattice_vectors
                        )
                        i0 = i1

        self.cell.lattice_vectors = np.array(
            [m[c] * self.cell.lattice_vectors[c] for c in range(self.dim)]
        )

        self.shape = tuple(self.shape * np.array(m))
        return self

    def draw(
        self,
        radius=None,
        show_labels=False,
        colouring=None,
        units=None,
        equal_aspect=True,
    ):
        """
        Visualize the positions of a set of qubits.

        Parameters
        ----------
        radius: float | str
            A cutoff radius for visualizing bonds.
            Pass 'nearest' to set the radius to the smallest
            distance between the passed qubits.
            If no value is passed the bonds will not be visualized.
        show_labels: bool
            Whether to show the labels of the qubits.
            Defaults to False.
        colouring: str | list
            A set of integers used to assign different colors to each Qubit.
            This can be used to view different magnetic orderings.
            Must have the same length as the number of Qubits.
        units : str, optional
            The units of distance.
        equal_aspect : bool, optional
            Whether to have the same scaling for the axes.
            Defaults to True.

        See Also
        --------
        qse.draw
        """
        _draw(
            self,
            radius=radius,
            show_labels=show_labels,
            colouring=colouring,
            units=units,
            equal_aspect=equal_aspect,
        )

    def repeat(self, rep):
        """Create new repeated qbits object.

        The *rep* argument should be a sequence of three positive
        integers like *(2,3,1)* or a single integer (*r*) equivalent
        to *(r,r,r)*."""

        qbits = self.copy()
        qbits *= rep
        return qbits

    def __mul__(self, rep):
        return self.repeat(rep)

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
        r"""
        Get the centroid of the positions.

        Returns
        -------
        np.ndarray
            The centroid of the positions.

        Notes
        -----
        For a set of :math:`k` positions
        :math:`\textbf{x}_1, \textbf{x}_2, ..., \textbf{x}_k`
        the centroid is given by

        .. math::

            \frac{\textbf{x}_1 + \textbf{x}_2 + ... + \textbf{x}_k}{k}.
        """
        return self.positions.mean(0)

    def set_centroid(self, centroid):
        r"""
        Set the centroid of the positions.

        Parameters
        ----------
        centroid : float | np.ndarray
            The new centroid. Can be a float or a xyz vector

        Notes
        -----
        For a set of :math:`k` positions
        :math:`\textbf{x}_1, \textbf{x}_2, ..., \textbf{x}_k`
        the centroid is given by

        .. math::

            \frac{\textbf{x}_1 + \textbf{x}_2 + ... + \textbf{x}_k}{k}.
        """
        self.positions += centroid - self.get_centroid()

    def rotate(self, a, v="z", center=(0, 0, 0), rotate_cell=False):
        r"""
        Rotate qbits based on a vector and an angle, or two vectors.

        Parameters
        ----------
        a :
            Angle that the qbits is rotated (anticlockwise) around the vector 'v'.
            'a' can also be a vector and then 'a' is rotated
            into 'v'. If 'a' is an angle it must be in degrees.
        v :
            Vector to rotate the qbits around. Vectors can be given as
            strings: 'x', '-x', 'y', ... .
            Defaults to 'z'.
        center :
            The center is kept fixed under the rotation. Use 'COP' to
            fix the center of positions or 'COU' to fix the center of
            cell. Defaults to = (0, 0, 0).
        rotate_cell = False:
            If true the cell is also rotated.

        Examples
        --------
        The following all rotate 90 degrees anticlockwise around the z-axis,
        so that the x-axis is rotated into the y-axis:

        >>> qbits.rotate(90)
        >>> qbits.rotate(90, 'z')
        >>> qbits.rotate(90, (0, 0, 1))
        >>> qbits.rotate(-90, '-z')
        >>> qbits.rotate('x', 'y')
        >>> qbits.rotate((1, 0, 0), (0, 1, 0))

        Notes
        -----
        If 'a' is an angle, :math:`\theta`, and if :math:`\textbf{v}` is the vector
        then we define

        .. math::
            R = \cos(\theta)I + \sin(\theta)[\textbf{v}]_\times
            + (1-\cos(\theta))\textbf{v}\textbf{v}^T

        where :math:`[\textbf{v}]_\times \textbf{x} = \textbf{v} \times \textbf{x}`.
        If
        :math:`\textbf{r}` is a coordinate vector
        and :math:`\textbf{c}` is the center, this transforms
        the coordinate vector to

        .. math::

            \textbf{r} \rightarrow R(\textbf{r}-\textbf{c}) + \textbf{c}.
        """
        if self.dim == 1:
            raise Exception("rotate requires 2 or 3 dimensions.")

        rm_dim = False
        if self.dim == 2:
            self.add_dim()
            rm_dim = True

        if not isinstance(a, numbers.Real):
            a, v = v, a

        v = _norm_vector(_string2vector(v))

        if isinstance(a, numbers.Real):
            a = _to_rads(a)
            c = cos(a)
            s = sin(a)
        else:
            v2 = _norm_vector(_string2vector(a))
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
        if rm_dim:
            self.remove_dim("z")

        if rotate_cell:
            self.cell.lattice_vectors = (
                c * self.cell.lattice_vectors
                - np.cross(self.cell.lattice_vectors, s * v)
                + np.outer(np.dot(self.cell.lattice_vectors, v), (1.0 - c) * v)
            )

    def _centering_as_array(self, center):
        if isinstance(center, str):
            if center.lower() == "cop":
                center = self.get_centroid()
            elif center.lower() == "cou":
                center = self.cell.lattice_vectors.sum(axis=0) / 2
            else:
                raise ValueError("Cannot interpret center")
        else:
            center = np.array(center, float)
        return center

    def euler_rotate(self, phi=0.0, theta=0.0, psi=0.0, center=(0, 0, 0)):
        r"""
        Rotate qbits via Euler angles (in degrees).

        Parameters
        ----------
        phi : float
            The 1st rotation angle around the z axis.
        theta : float
            Rotation around the x axis.
        psi : float
            2nd rotation around the z axis.
        center :
            The point to rotate about. A sequence of length 3 with the
            coordinates, or 'COM' to select the center of mass, 'COP' to
            select center of positions or 'COU' to select center of cell.

        Notes
        -----
        Let

        .. math::

            R =
            \begin{pmatrix}
            \cos(\psi ) & \sin(\psi ) & 0 \\
            -\sin(\psi ) & \cos(\psi ) & 0\\
            0 & 0 & 1\\
            \end{pmatrix}
            \begin{pmatrix}
            1 & 0 & 0\\
            0 & \cos(\theta ) & \sin(\theta ) \\
            0 & -\sin(\theta ) & \cos(\theta ) \\
            \end{pmatrix}
            \begin{pmatrix}
            \cos(\phi ) & \sin(\phi ) & 0 \\
            -\sin(\phi ) & \cos(\phi ) & 0\\
            0 & 0 & 1\\
            \end{pmatrix}

        then if :math:`\textbf{r}` is a coordinate vector
        and :math:`\textbf{c}` is the center, this transforms
        the coordinate vector to

        .. math::

            \textbf{r} \rightarrow R(\textbf{r}-\textbf{c}) + \textbf{c}.
        """
        if self.dim != 3:
            raise Exception("euler_rotate can only be performed on 3D systems.")

        def rotation_mat(angle):
            return np.array(
                [
                    [np.cos(angle), np.sin(angle)],
                    [-np.sin(angle), np.cos(angle)],
                ]
            )

        # First Euler rotation about z in matrix form
        D = np.eye(3)
        D[:-1, :-1] = rotation_mat(_to_rads(phi))

        # Second Euler rotation about x:
        C = np.eye(3)
        C[1:, 1:] = rotation_mat(_to_rads(theta))

        # Third Euler rotation, 2nd rotation about z:
        B = np.eye(3)
        B[:-1, :-1] = rotation_mat(_to_rads(psi))

        # Total Euler rotation
        A = np.dot(B, np.dot(C, D))

        # Move the molecule to the origin.
        rcoords = self.positions - self._centering_as_array(center)

        # Do the rotation
        rcoords = np.dot(A, np.transpose(rcoords))

        # Move back to the rotation point
        self.positions = np.transpose(rcoords) + center

    def _masked_rotate(self, center, axis, diff, mask):
        # do rotation of subgroup by copying it to temporary qbits object
        # and then rotating that
        #
        # recursive object definition might not be the most elegant thing,
        # more generally useful might be a rotation function with a mask?
        group = self.__class__()
        for i in range(len(self)):
            if mask[i]:
                group += self[i]
        group.translate(-center)
        group.rotate(diff * 180 / np.pi, axis)
        group.translate(center)
        # set positions in original qbits object
        j = 0
        for i in range(len(self)):
            if mask[i]:
                self.positions[i] = group[j].position
                j += 1

    def get_angle(self, i: int, j: int, k: int):
        """
        Get the angle in degress formed by three qubits.

        Parameters
        ----------
        i : int
            The index of the first qubit.
        j : int
            The index of the second qubit.
        k : int
            The index of the third qubit.

        Returns
        -------
        float
            The angle between the qubits.

        Notes
        -----
        Let x1, x2, x3 be the vectors describing the positions of the three
        qubits. Then we calcule the angle between x1-x2 and x3-x2.
        """
        v1 = _norm_vector(self.positions[i] - self.positions[j])
        v2 = _norm_vector(self.positions[k] - self.positions[j])
        dot_prod = np.dot(v1, v2)

        # The if-statements are in case of floating point errors.
        if dot_prod >= 1.0:
            return 0.0
        if dot_prod <= -1.0:
            return 180.0
        return _to_degrees(np.arccos(dot_prod))

    def get_angles(self, indices):
        """
        Get the angle in degress formed by three qubits for multiple groupings.

        Parameters
        ----------
        indices : list | np.ndarray
            The indices of the groupings of qubits.
            Must be of shape (n, 3), where n is the number of groupings.

        Returns
        -------
        np.ndarray
            The angles between the qubits.

        Notes
        -----
        Let x1, x2, x3 be the vectors describing the positions of the three
        qubits. Then we calcule the angle between x1-x2 and x3-x2 for all the
        different groupings.
        """
        indices = np.array(indices)
        if indices.shape[1] != 3:
            raise Exception("The indicies must be of shape (-1, 3).")
        return np.array([self.get_angle(i, j, k) for i, j, k in indices])

    def set_angle(
        self, a1, a2=None, a3=None, angle=None, mask=None, indices=None, add=False
    ):
        """
        Set angle (in degrees) formed by three qbits.

        Sets the angle between vectors *a2*->*a1* and *a2*->*a3*.

        If *add* is `True`, the angle will be changed by the value given.

        Same usage as in :meth:`ase.Qbits.set_dihedral`.
        If *mask* and *indices*
        are given, *indices* overwrites *mask*. If *mask* and *indices*
        are not set, only *a3* is moved.
        """

        if any(a is None for a in [a2, a3, angle]):
            raise ValueError("a2, a3, and angle must not be None")

        # If not provided, set mask to the last qbit in the angle description
        if mask is None and indices is None:
            mask = np.zeros(len(self))
            mask[a3] = 1
        elif indices is not None:
            mask = [index in indices for index in range(len(self))]

        if add:
            diff = angle
        else:
            # Compute necessary in angle change, from current value
            diff = angle - self.get_angle(a1, a2, a3)

        diff = _to_rads(diff)
        # Do rotation of subgroup by copying it to temporary qbits object and
        # then rotating that
        v10 = self.positions[a1] - self.positions[a2]
        v12 = self.positions[a3] - self.positions[a2]
        v10 /= np.linalg.norm(v10)
        v12 /= np.linalg.norm(v12)
        axis = np.cross(v10, v12)
        center = self.positions[a2]
        self._masked_rotate(center, axis, diff, mask)

    def rattle(self, stdev=0.001, seed=None, rng=None):
        """
        Randomly displace qbits.

        This method adds random displacements to the qbit positions.
        The random numbers are
        drawn from a normal distribution of standard deviation stdev.

        For a parallel calculation, it is important to use the same
        seed on all processors!
        """

        if seed is not None and rng is not None:
            raise ValueError("Please do not provide both seed and rng.")

        if rng is None:
            if seed is None:
                seed = 42
            rng = np.random.RandomState(seed)
        self.positions += rng.normal(scale=stdev, size=self.positions.shape)

    def get_distance(self, i, j):
        """
        Return the distance between two qbits.

        Parameters
        ----------
        i : int
            The index of the first qubit.
        j : int
            The index of the second qubit.

        Returns
        -------
        float
            The distance between the qubits.
        """
        return np.linalg.norm(self.positions[i] - self.positions[j])

    def get_distances(self, i, indices):
        """
        Return distances of the ith qubit with a list of qubits.

        Parameters
        ----------
        i : int
            The index of the ith qubit.
        indices : list[int]
            The indices of other qubits.

        Returns
        -------
        np.ndarray
            An array containing the distances.
        """
        return np.array([self.get_distance(i, j) for j in indices])

    def get_all_distances(self):
        """
        Return the distances of all of the qubits with all of the other qubits.

        Returns
        -------
        np.ndarray
            An array of shape (nqbits, nqbits) containing the distances.
        """
        distances = np.zeros((self.nqbits, self.nqbits))
        for i in range(self.nqbits - 1):
            for j in range(i + 1, self.nqbits):
                distances[i, j] = distances[j, i] = self.get_distance(i, j)

        return distances

    def set_distance(
        self,
        i,
        j,
        distance,
        fix=0.5,
        mask=None,
        indices=None,
        add=False,
        factor=False,
    ):
        """
        Set the distance between qubits i and j.

        Parameters
        ----------
        i : int
            The index of the ith qubit.
        j : int
            The index of the jth qubit.
        distance : float
            The new distance to be set.
        fix : float
            By default, the center of the two qbits will be fixed.  Use
            fix=0 to fix the first qbit, fix=1 to fix the second
            qbit and fix=0.5 (default) to fix the center of the bond.
        mask: np.ndarray | list
            If mask or indices are set (mask overwrites indices),
            only the qbits defined there are moved. It is assumed that the
            qbits in mask/indices move together with the jth qubit.
            If fix=1, only the ith qubit will therefore be moved.
        indices: np.ndarray | list
            If mask or indices are set (mask overwrites indices),
            only the qbits defined there are moved. It is assumed that the
            qbits in mask/indices move together with the jth qubit.
            If fix=1, only the ith qubit will therefore be moved.
        add:
            When add is true, the distance is changed by the value given.
        factor:
            When factor is true, the value given is a factor scaling the distance.
        """

        if i % len(self) == j % len(self):
            raise ValueError("i and j must not be the same")

        if add:
            distance += self.get_distance(i, j)
        elif factor:
            distance *= self.get_distance(i, j)

        if mask is None and indices is None:
            indices = [i, j]
        elif mask:
            indices = [ind for ind in range(len(self)) if mask[ind]]

        distance_vec = self.positions[j] - self.positions[i]
        x = 1.0 - distance / np.linalg.norm(distance_vec)

        for ind in indices:
            if ind == i:
                self.positions[ind] += (x * fix) * distance_vec
            else:
                self.positions[ind] -= (x * (1.0 - fix)) * distance_vec

    def wrap(self, **wrap_kw):
        """Wrap positions to unit cell.

        Parameters:

        wrap_kw: (keyword=value) pairs
            optional keywords `pbc`, `center`, `pretty_translation`, `eps`,
            see :func:`ase.geometry.wrap_positions`
        """

        if "pbc" not in wrap_kw:
            wrap_kw["pbc"] = self.pbc

        self.positions[:] = self.get_positions(wrap=True, **wrap_kw)

    def __eq__(self, other):
        """Check for identity of two qbits objects.

        Identity means: same positions, states, unit cell and
        periodic boundary conditions."""
        if not isinstance(other, Qbits):
            return False
        a = self.arrays
        b = other.arrays
        return (
            len(self) == len(other)
            and (a["positions"] == b["positions"]).all()
            and (a["states"] == b["states"]).all()
            and (self.cell == other.cell).all()
            and (self.pbc == other.pbc).all()
        )

    def __ne__(self, other):
        """Check if two qbits objects are not equal.

        Any differences in positions, states, unit cell or
        periodic boundary condtions make qbits objects not equal.
        """
        eq = self.__eq__(other)
        if eq is NotImplemented:
            return eq
        else:
            return not eq

    # Rajarshi: Below these three written to add attribute of states
    def _get_states(self):
        """Return reference to states-array for in-place manipulations."""
        return self.arrays["states"]

    def _set_states(self, sts):
        """Set states directly, bypassing constraints."""
        self.arrays["states"][
            :
        ] = sts  # (sts.T / np.linalg.norm(sts, axis=1)).T # need to be normalized

    # below is equivalent to defining @property states
    states = property(
        _get_states,
        _set_states,
        doc="Attribute for direct " + "manipulation of the states.",
    )

    # Rajarshi: Write method to get labels
    def _get_labels(self):
        """Return array of labels"""
        return self.arrays["labels"]

    def _set_labels(self, lbs):
        """Set the labels directly."""
        self.arrays["labels"][:] = lbs

    labels = property(
        _get_labels, _set_labels, doc="Attribute for direct manipulation of labels"
    )

    def write(self, filename, format=None, **kwargs):
        """Write qbits object to a file.

        see ase.io.write for formats.
        kwargs are passed to ase.io.write.
        """
        from pickle import dump

        dump(obj=self.todict(), file=open(filename, "wb"), **kwargs)

    def iterimages(self):
        yield self

    def to_pulser(self):
        from pulser import Register

        if self.dim == 1:
            warnings.warn("1D system passed, adding a y axis.")
            return Register.from_coordinates(
                np.column_stack([self.positions, np.zeros(self.nqbits)]), prefix="q"
            )

        if self.dim == 2:
            return Register.from_coordinates(self.positions, prefix="q")

        if self.dim == 3:
            warnings.warn("3D system passed, removing the z axis.")
            return Register.from_coordinates(self.positions[:, :2], prefix="q")

        return Exception("The qbits must be 2D or 3D for use in Pulser.")

    def compute_interaction_hamiltonian(
        self,
        distance_func,
        interaction,
        tol=1e-8,
    ):
        """
        Compute the interaction Hamiltonian for a system of qubits.

        This function constructs an Operators object based on the distances between
        the qubits.

        Parameters
        ----------
        distance_func : callable
            A function that takes a distance (float) and returns the interaction
            coefficient (float).
        interaction : str | list[str]
            The type of interaction (e.g., "X", "Y", "Z") for the Hamiltonian terms,
            or can be a list of strings.
        tol : float, optional
            Tolerance threshold for including interaction terms. Terms with absolute
            coefficients less than `tol` are discarded. Default is 1e-8.

        Returns
        -------
        Operators
            The interaction operators.

        Examples
        --------
        To create a ZZ Hamiltonian for only nearest neighbour qubits

        >>> spacing = 1.0
        >>> qbits = qse.lattices.chain(spacing, 4)
        >>> coupling = -2.
        >>> qbits.compute_interaction_hamiltonian(
        ...     lambda x: coupling*np.isclose(x, spacing), "Z"
        ... )
        ... Number of qubits: 4
        ... Number of terms: 3
        ...
        ... -2.00 Z0 Z1
        ... -2.00 Z1 Z2
        ... -2.00 Z2 Z3

        To create an XY Hamiltonian based on distance

        >>> spacing = 1.0
        >>> qbits = qse.lattices.chain(spacing, 2)
        >>> coupling = 1.
        >>> hamiltonian = qbits.compute_interaction_hamiltonian(
        ...     lambda x: coupling / x**3, ["X", "Y"]
        ... )
        >>> hamiltonian += qbits.compute_interaction_hamiltonian(
        ...     lambda x: coupling / x**3, ["Y", "X"]
        ... )
        ... Number of qubits: 2
        ... Number of terms: 2
        ...
        ... 1.00 X0 Y1
        ... 1.00 Y0 X1
        """
        ops = []

        for i in range(self.nqbits - 1):
            for j in range(i + 1, self.nqbits):
                coef = distance_func(self.get_distance(i, j))
                if np.abs(coef) > tol:
                    ops.append(Operator(interaction, (i, j), self.nqbits, coef))

        return Operators(ops)

    def add_dim(self):
        """
        Adds a spatial dimension to the positions. E.e. to go from 1D to 2D systems
        or 2D to 3D systems.
        """
        if self.dim == 3:
            raise ValueError("Can't go above 3 dimensions.")
        self.arrays["positions"] = np.column_stack(
            [self.arrays["positions"], np.zeros(self.nqbits)]
        )

    def remove_dim(self, dim):
        """
        Removes a spatial dimension to the positions. E.e. to go from 2D to 1D systems
        or 3D to 2D systems.

        Parameters
        ----------
        dim : str
            The dimension to be removed.
            Must be one of 'x', 'y' or 'z'.
        """
        if self.dim == 1:
            raise ValueError("Can't go below 1 dimension.")
        axes = ["x", "y", "z"]
        if dim not in axes:
            raise ValueError("dim must be one of 'x', 'y' or 'z'.")
        keep_cols = [i for i in range(self.dim) if i != axes.index(dim)]

        self.arrays["positions"] = self.arrays["positions"][:, keep_cols]

    def get_scaled_positions(self):
        """
        Get the positions in units of the unit cell.

        Returns
        -------
        np.ndarray
            The positions in units of the unit cell.
        """
        return np.dot(self.positions, np.linalg.inv(self.cell.lattice_vectors))


def _norm_vector(v):
    normv = np.linalg.norm(v)

    if normv == 0.0:
        raise ZeroDivisionError("Vector has 0 norm.", v)
    return v / normv


def _to_rads(a):
    return a * np.pi / 180


def _to_degrees(a):
    return 180 * a / np.pi


def _string2vector(v):
    """Used in rotate method to rotate qbit location"""
    if isinstance(v, str):
        if v[0] == "-":
            return -_string2vector(v[1:])
        w = np.zeros(3)
        w["xyz".index(v)] = 1.0
        return w
    return np.array(v, float)


def default(data, dflt):
    """Helper function for setting default values."""
    if data is None:
        return None
    elif isinstance(data, (list, tuple)):
        newdata = []
        allnone = True
        for x in data:
            if x is None:
                newdata.append(dflt)
            else:
                newdata.append(x)
                allnone = False
        if allnone:
            return None
        return newdata
    else:
        return data
