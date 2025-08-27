"""
This module defines the central object in the QSE package: the Qbits object.
"""

import numbers
from math import cos, sin

import numpy as np
from ase.cell import Cell
from ase.geometry import (
    find_mic,
)

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
    labels: list of str
        A list of strings corresponding to a label for each qubit.
    states: list of 2-length arrays.
        State of each qubit.
    positions: list of xyz-positions
        Qubit positions.  Anything that can be converted to an
        ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
        ...].
    scaled_positions: list of scaled-positions
        Like positions, but given in units of the unit cell.
        Can not be set at the same time as positions.
    cell: 3x3 matrix or length 3 or 6 vector
        Unit cell vectors.  Can also be given as just three
        numbers for orthorhombic cells, or 6 numbers, where
        first three are lengths of unit cell vectors, and the
        other three are angles between them (in degrees), in following order:
        [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)].
        First vector will lie in x-direction, second in xy-plane,
        and the third one in z-positive subspace.
        Default value: [0, 0, 0].
    calculator: calculator object
        Used to attach a calculator for doing computation.

    Examples
    --------
    Empty Qbits object:

    >>> qs = qse.Qbits()

    These are equivalent:

    >>> a = qse.Qbits(
    ...     labels=['qb1', 'qb2'],
    ...     positions=np.array([(0, 0, 0), (0, 0, 2)])
    ... )
    >>> a = qse.Qbits.from_qbit_list(
    ...     [Qbit('qb1', position=(0, 0, 0)), Qbit('qb2', position=(0, 0, 2))]
    ... )

    >>> xd = np.array(
    ...    [[0, 0, 0],
    ...     [0.5, 0.5, 0.5]])
    >>> qdim = qse.Qbits(positions=xd)
    >>> qdim.cell = [1,1,1]
    >>> qlat = qdim.repeat([3,3,3])

    The qdim will have shape = (2,1,1) and qlat will have shape = (6, 3, 3)

    Notes
    -----
    In order to do computation, a calculator object has to attached
    to the qbits object.
    """

    def __init__(
        self,
        labels=None,
        states=None,
        positions=None,
        scaled_positions=None,
        cell=None,
        calculator=None,
    ):
        if (positions is not None) and (scaled_positions is not None):
            raise Exception(
                "Both 'positions' and 'scaled_positions'"
                " cannot be passed at the same time."
            )

        if (scaled_positions is not None) and (cell is None):
            raise Exception("'scaled_positions' requires 'cell' to not be None.")

        # get number of qubits
        if labels is None:
            if positions is not None:
                nqbits = len(positions)
            elif scaled_positions is not None:
                nqbits = len(scaled_positions)
            else:
                nqbits = 0
        else:
            if not isinstance(labels, list):
                raise Exception("'labels' must be a list.")
            nqbits = len(labels)

        if (positions is not None) and (len(positions) != nqbits):
            raise Exception("Both 'positions' and 'labels' must have the same length.")

        if (scaled_positions is not None) and (len(scaled_positions) != nqbits):
            raise Exception(
                "Both 'scaled_positions' and 'labels' must have the same length."
            )

        self.arrays = {}

        # labels
        if labels is None:
            labels = [str(i) for i in range(nqbits)]
        # We allow for labels up to length 12.
        self.new_array("labels", labels, "<U12")

        # cell
        self._cellobj = Cell.new()
        if cell is None:
            cell = np.zeros((3, 3))
        self.set_cell(cell)

        # positions
        if positions is None:
            if scaled_positions is None:
                positions = np.zeros((len(self.arrays["labels"]), 3))
            else:
                assert self.cell.rank == 3
                positions = np.dot(scaled_positions, self.cell)
        self.new_array("positions", positions, float, (3,))

        # states
        if states is None:
            states = np.zeros((positions.shape[0], 2), dtype=complex)
            states[:, 0] = 1
        self.new_array("states", states, complex, (2,))

        # shape
        self._shape = (self.nqbits, 1, 1)

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
    def shape(self):
        """The shape of the qbits"""
        return self._shape

    @shape.setter
    def shape(self, new_shape):
        """Update the shape to new shape"""
        if self.nqbits != np.prod(new_shape):
            raise AssertionError(f"no. of qubits= {self.nqbits}, yet shape {new_shape}")
        self._shape = new_shape

    @property
    def calc(self):
        """Calculator object."""
        return self._calc

    @calc.setter
    def calc(self, calc):
        self._calc = calc
        if hasattr(calc, "set_qbits"):
            calc.set_qbits(self)

    def set_cell(self, cell, scale_qbits=False):
        """Set unit cell vectors.

        Parameters:

        cell: 3x3 matrix or length 3 or 6 vector
            Unit cell.  A 3x3 matrix (the three unit cell vectors) or
            just three numbers for an orthorhombic cell. Another option is
            6 numbers, which describes unit cell with lengths of unit cell
            vectors and with angles between them (in degrees), in following
            order: [len(a), len(b), len(c), angle(b,c), angle(a,c),
            angle(a,b)].  First vector will lie in x-direction, second in
            xy-plane, and the third one in z-positive subspace.
        scale_qbits: bool
            Fix qbit positions or move qbits with the unit cell?
            Default behavior is to *not* move the qbits (scale_qbits=False).

        Examples:

        Two equivalent ways to define an orthorhombic cell:

        >>> qbits = Qbits('He')
        >>> a, b, c = 7, 7.5, 8
        >>> qbits.set_cell([a, b, c])
        >>> qbits.set_cell([(a, 0, 0), (0, b, 0), (0, 0, c)])

        FCC unit cell:

        >>> qbits.set_cell([(0, b, b), (b, 0, b), (b, b, 0)])

        Hexagonal unit cell:

        >>> qbits.set_cell([a, a, c, 90, 90, 120])

        Rhombohedral unit cell:

        >>> alpha = 77
        >>> qbits.set_cell([a, a, a, alpha, alpha, alpha])
        """

        # print('here', cell)
        cell = Cell.new(cell)

        if scale_qbits:
            M = np.linalg.solve(self.cell.complete(), cell.complete())
            self.positions[:] = np.dot(self.positions, M)

        self.cell[:] = cell

    def get_cell(self, complete=False):
        """Get the three unit cell vectors as a `class`:ase.cell.Cell` object.

        The Cell object resembles a 3x3 ndarray, and cell[i, j]
        is the jth Cartesian coordinate of the ith cell vector."""
        if complete:
            return self.cell.complete()
        return self.cell.copy()

    @property
    def nqbits(self):
        return len(self)

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
        """Update array.

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

    def has(self, name):
        """
        Check for existence of array.

        name must be one of: 'momenta', 'masses', 'initial_magmoms',
        'initial_charges'.
        """
        # XXX extend has to calculator properties
        return name in self.arrays

    def copy(self):
        """Return a copy."""
        qbits = self.__class__(cell=self.cell)

        qbits.arrays = {}
        for name, a in self.arrays.items():
            qbits.arrays[name] = a.copy()

        qbits.shape = self.shape  # this was necessary, and took long time to realise!

        return qbits

    def todict(self):
        """For basic JSON (non-database) support."""
        # d = dict(self.arrays)
        d = {}
        d["labels"] = self.arrays["labels"]
        d["positions"] = self.arrays["positions"]
        d["states"] = self.arrays["states"]
        d["cell"] = self.cell  # np.asarray(self.cell)
        # Calculator...  trouble.
        return d

    @classmethod
    def fromdict(cls, dct):
        """Rebuild qbits object from dictionary representation (todict)."""
        dct = dct.copy()
        kw = {}
        for name in ["labels", "positions", "states", "cell"]:
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

        cell = self.cell
        if cell:
            if cell.orthorhombic:
                cell = cell.lengths().tolist()
            else:
                cell = cell.tolist()
            tokens.append("cell={0}".format(cell) + ",\n")
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

        n1 = len(self)
        n2 = len(other)

        for name, a1 in self.arrays.items():
            a = np.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            if name == "masses":
                a2 = other.get_masses()
            else:
                a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = np.empty((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            if name == "masses":
                a[:n1] = self.get_masses()[:n1]
            else:
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
            is a list or a slice, a Qbits object with the same cell, and
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

        qbits = self.__class__(
            cell=self.cell,
        )

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
            m = (m, m, m)

        for x, vec in zip(m, self.cell):
            if x != 1 and not vec.any():
                raise ValueError("Cannot repeat along undefined lattice " "vector")

        M = np.prod(m)
        n = len(self)

        for name, a in self.arrays.items():
            self.arrays[name] = np.tile(a, (M,) + (1,) * (len(a.shape) - 1))

        positions = self.arrays["positions"]
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), self.cell)
                    i0 = i1

        self.cell = np.array([m[c] * self.cell[c] for c in range(3)])

        new_shape = tuple(self.shape * np.array(m))
        self.shape = new_shape
        return self

    def draw(self, radius=None, show_labels=False, colouring=None, units=None):
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
        if rotate_cell:
            rotcell = self.get_cell()
            rotcell[:] = (
                c * rotcell
                - np.cross(rotcell, s * v)
                + np.outer(np.dot(rotcell, v), (1.0 - c) * v)
            )
            self.set_cell(rotcell)

    def _centering_as_array(self, center):
        if isinstance(center, str):
            if center.lower() == "cop":
                center = self.get_centroid()
            elif center.lower() == "cou":
                center = self.get_cell().sum(axis=0) / 2
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
        """Randomly displace qbits.

        This method adds random displacements to the qbit positions.
        The random numbers are drawn from a normal distribution of
        standard deviation stdev.

        For a parallel calculation, it is important to use the same
        seed on all processors!"""

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
        a0,
        a1,
        distance,
        fix=0.5,
        mic=False,
        mask=None,
        indices=None,
        add=False,
        factor=False,
    ):
        """Set the distance between two qbits.

        Set the distance between qbits *a0* and *a1* to *distance*.
        By default, the center of the two qbits will be fixed.  Use
        *fix=0* to fix the first qbit, *fix=1* to fix the second
        qbit and *fix=0.5* (default) to fix the center of the bond.

        If *mask* or *indices* are set (*mask* overwrites *indices*),
        only the qbits defined there are moved
        (see :meth:`ase.Qbits.set_dihedral`).

        When *add* is true, the distance is changed by the value given.
        In combination
        with *factor* True, the value given is a factor scaling the distance.

        It is assumed that the qbits in *mask*/*indices* move together
        with *a1*. If *fix=1*, only *a0* will therefore be moved."""

        if a0 % len(self) == a1 % len(self):
            raise ValueError("a0 and a1 must not be the same")

        if add:
            oldDist = self.get_distance(a0, a1, mic=mic)
            if factor:
                newDist = oldDist * distance
            else:
                newDist = oldDist + distance
            self.set_distance(
                a0,
                a1,
                newDist,
                fix=fix,
                mic=mic,
                mask=mask,
                indices=indices,
                add=False,
                factor=False,
            )
            return

        R = self.arrays["positions"]
        D = np.array([R[a1] - R[a0]])

        if mic:
            D, D_len = find_mic(D, self.cell)
        else:
            D_len = np.array([np.sqrt((D**2).sum())])
        x = 1.0 - distance / D_len[0]

        if mask is None and indices is None:
            indices = [a0, a1]
        elif mask:
            indices = [i for i in range(len(self)) if mask[i]]

        for i in indices:
            if i == a0:
                R[a0] += (x * fix) * D[0]
            else:
                R[i] -= (x * (1.0 - fix)) * D[0]

    def set_scaled_positions(self, scaled):
        """Set positions relative to unit cell."""
        self.positions[:] = self.cell.cartesian_positions(scaled)

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

    def get_volume(self):
        """Get volume of unit cell."""
        if self.cell.rank != 3:
            raise ValueError(
                "You have {0} lattice vectors: volume not defined".format(
                    self.cell.rank
                )
            )
        return self.cell.volume

    def _get_positions(self):
        """Return reference to positions-array for in-place manipulations."""
        return self.arrays["positions"]

    def _set_positions(self, pos):
        """Set positions directly."""
        self.arrays["positions"][:] = pos

    positions = property(
        _get_positions,
        _set_positions,
        doc="Attribute for direct " + "manipulation of the positions.",
    )

    # Rajarshi: Below these three written to add attribute of states
    def _get_states(self):
        """Return reference to states-array for in-place manipulations."""
        return self.arrays["states"]

    def _set_states(self, sts):
        """Set states directly."""
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

    @property
    def cell(self):
        """The :class:`ase.cell.Cell` for direct manipulation."""
        return self._cellobj

    @cell.setter
    def cell(self, cell):
        cell = Cell.ascell(cell)
        self._cellobj[:] = cell

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

        return Register.from_coordinates(self.positions[:, :2], prefix="q")

    #

    # Rajarshi: Deleted the edit method, which in original
    # ASE approach lets users manipulate Atoms object. At
    # some stage we may adopt similar approach depending on
    # the usage/usecase.
    # def edit(self): Modify qbits interactively through ASE's GUI viewer.


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
