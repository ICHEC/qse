# Copyright 2008, 2009 CAMd
# (see accompanying license files for details).

"""
Definition of the Qbits class.

This module defines the central object in the QSE package: the Qbits object.
"""
import copy
import numbers

import numpy as np
from ase.cell import Cell
from ase.geometry import (
    find_mic,
    get_angles,
    get_dihedrals,
    get_distances,
    wrap_positions,
)
from ase.utils import deprecated

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
    celldisp: Vector
        Unit cell displacement vector. To visualize a displaced cell
        around the center of mass of a Systems of qbits. Default value
        = (0,0,0)
    pbc: one or three bool
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        value: False.
    constraint: constraint object(s)
        Used for applying one or more constraints during structure
        optimization.
    calculator: calculator object
        Used to attach a calculator for doing computation.
    info: dict of key-value pairs
        Dictionary of key-value pairs with additional information
        about the system.  The following keys may be used by ase:

          - spacegroup: Spacegroup instance
          - unit_cell: 'conventional' | 'primitive' | int | 3 ints
          - adsorbate_info: Information about special adsorption sites

        Items in the info attribute survives copy and slicing and can
        be stored in and retrieved from trajectory files given that the
        key is a string, the value is JSON-compatible and, if the value is a
        user-defined object, its base class is importable.  One should
        not make any assumptions about the existence of keys.

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
    >>> qdim.pbc = True
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
        pbc=None,
        celldisp=None,
        constraint=None,
        calculator=None,
        info=None,
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

        # pbc
        self._pbc = np.zeros(3, bool)
        if pbc is None:
            pbc = False
        self.set_pbc(pbc)

        # celldisp
        if celldisp is None:
            celldisp = np.zeros(shape=(3, 1))
        self.set_celldisp(celldisp)

        # constraint
        self.set_constraint(constraint)

        # calculator
        self.calc = calculator

        # info
        self.info = {} if info is None else dict(info)

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

    def set_constraint(self, constraint=None):
        """Apply one or more constrains.

        The *constraint* argument must be one constraint object or a
        list of constraint objects."""
        if constraint is None:
            self._constraints = []
        else:
            if isinstance(constraint, list):
                self._constraints = constraint
            elif isinstance(constraint, tuple):
                self._constraints = list(constraint)
            else:
                self._constraints = [constraint]

    def _get_constraints(self):
        return self._constraints

    def _del_constraints(self):
        self._constraints = []

    constraints = property(
        _get_constraints, set_constraint, _del_constraints, "Constraints of the qbits."
    )

    def set_cell(self, cell, scale_qbits=False, apply_constraint=True):
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
        apply_constraint: bool
            Whether to apply constraints to the given cell.

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

        # Override pbcs if and only if given a Cell object:
        # print('here', cell)
        cell = Cell.new(cell)

        # XXX not working well during initialize due to missing _constraints
        if apply_constraint and hasattr(self, "_constraints"):
            for constraint in self.constraints:
                if hasattr(constraint, "adjust_cell"):
                    constraint.adjust_cell(self, cell)

        if scale_qbits:
            M = np.linalg.solve(self.cell.complete(), cell.complete())
            self.positions[:] = np.dot(self.positions, M)

        self.cell[:] = cell

    def set_celldisp(self, celldisp):
        """Set the unit cell displacement vectors."""
        celldisp = np.array(celldisp, float)
        self._celldisp = celldisp

    def get_celldisp(self):
        """Get the unit cell displacement vectors."""
        return self._celldisp.copy()

    def get_cell(self, complete=False):
        """Get the three unit cell vectors as a `class`:ase.cell.Cell` object.

        The Cell object resembles a 3x3 ndarray, and cell[i, j]
        is the jth Cartesian coordinate of the ith cell vector."""
        if complete:
            return self.cell.complete()
        return self.cell.copy()

    @deprecated("Please use qbits.cell.cellpar() instead")
    def get_cell_lengths_and_angles(self):
        """Get unit cell parameters. Sequence of 6 numbers.

        First three are unit cell vector lengths and second three
        are angles between them::

            [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]

        in degrees.
        """
        return self.cell.cellpar()

    @deprecated("Please use qbits.cell.reciprocal()")
    def get_reciprocal_cell(self):
        """Get the three reciprocal lattice vectors as a 3x3 ndarray.

        Note that the commonly used factor of 2 pi for Fourier
        transforms is not included here."""

        return self.cell.reciprocal()

    @property
    def nqbits(self):
        return len(self)

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

    def set_positions(self, newpositions, apply_constraint=True):
        """
        Set positions, honoring any constraints. To ignore constraints,
        use *apply_constraint=False*.
        """
        if self.constraints and apply_constraint:
            newpositions = np.array(newpositions, float)
            for constraint in self.constraints:
                constraint.adjust_positions(self, newpositions)

        self.set_array("positions", newpositions, shape=(3,))

    def get_positions(self, wrap=False, **wrap_kw):
        """
        Get array of positions.

        Parameters:

        wrap: bool
            wrap qbits back to the cell before returning positions
        wrap_kw: (keyword=value) pairs
            optional keywords `pbc`, `center`, `pretty_translation`, `eps`,
            see :func:`ase.geometry.wrap_positions`
        """
        if wrap:
            if "pbc" not in wrap_kw:
                wrap_kw["pbc"] = self.pbc
            return wrap_positions(self.positions, self.cell, **wrap_kw)
        else:
            return self.arrays["positions"].copy()

    def get_properties(self, properties):
        """This method is experimental; currently for internal use."""
        # XXX Something about constraints.
        if self._calc is None:
            raise RuntimeError("Qbits object has no calculator.")
        return self._calc.calculate_properties(self, properties)

    def copy(self):
        """Return a copy."""
        qbits = self.__class__(
            cell=self.cell, pbc=self.pbc, info=self.info, celldisp=self._celldisp.copy()
        )

        qbits.arrays = {}
        for name, a in self.arrays.items():
            qbits.arrays[name] = a.copy()
        qbits.constraints = copy.deepcopy(self.constraints)
        #
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
        d["pbc"] = self.pbc
        if self._celldisp.any():
            d["celldisp"] = self._celldisp
        if self.constraints:
            d["constraints"] = self.constraints
        if self.info:
            d["info"] = self.info
        # Calculator...  trouble.
        return d

    @classmethod
    def fromdict(cls, dct):
        """Rebuild qbits object from dictionary representation (todict)."""
        dct = dct.copy()
        kw = {}
        for name in ["labels", "positions", "states", "cell", "pbc"]:
            kw[name] = dct.pop(name)

        constraints = dct.pop("constraints", None)
        if constraints:
            from ase.constraints import dict2constraint

            constraints = [dict2constraint(d) for d in constraints]

        # labels = dct.pop('labels', None)

        info = dct.pop("info", None)

        qbits = cls(
            constraint=constraints, celldisp=dct.pop("celldisp", None), info=info, **kw
        )
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

        qbits = self.__class__(
            cell=self.cell,
            pbc=self.pbc,
            info=self.info,
            celldisp=self._celldisp,
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

        if self.constraints is not None:
            self.constraints = [c.repeat(m, n) for c in self.constraints]

        self.cell = np.array([m[c] * self.cell[c] for c in range(3)])

        new_shape = tuple(self.shape * np.array(m))
        self.shape = new_shape
        return self

    def draw(self, ax=None, radius=None):
        _draw(self, ax=ax, radius=radius)

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
        self.arrays["positions"] += np.array(displacement)

    def center_in_unit_cell(self, vacuum=None, axis=(0, 1, 2), about=None):
        """
        Center qbits in unit cell.

        Centers the qbits in the unit cell, so there is the same
        amount of vacuum on all sides.

        vacuum: float (default: None)
            If specified adjust the amount of vacuum when centering.
            If vacuum=10.0 there will thus be 10 Angstrom of vacuum
            on each side.
        axis: int or sequence of ints
            Axis or axes to act on.  Default: Act on all axes.
        about: float or array (default: None)
            If specified, center the qbits about <about>.
            I.e., about=(0., 0., 0.) (or just "about=0.", interpreted
            identically), to center about the origin.
        """

        # Find the orientations of the faces of the unit cell
        cell = self.cell.complete()
        dirs = np.zeros_like(cell)

        lengths = cell.lengths()
        for i in range(3):
            dirs[i] = np.cross(cell[i - 1], cell[i - 2])
            dirs[i] /= np.linalg.norm(dirs[i])
            if dirs[i] @ cell[i] < 0.0:
                dirs[i] *= -1

        if isinstance(axis, int):
            axes = (axis,)
        else:
            axes = axis

        # Now, decide how much each basis vector should be made longer
        pos = self.positions
        longer = np.zeros(3)
        shift = np.zeros(3)
        for i in axes:
            if len(pos):
                scalarprod = pos @ dirs[i]
                p0 = scalarprod.min()
                p1 = scalarprod.max()
            else:
                p0 = 0
                p1 = 0
            height = cell[i] @ dirs[i]
            if vacuum is not None:
                lng = (p1 - p0 + 2 * vacuum) - height
            else:
                lng = 0.0  # Do not change unit cell size!
            top = lng + height - p1
            shf = 0.5 * (top - p0)
            cosphi = cell[i] @ dirs[i] / lengths[i]
            longer[i] = lng / cosphi
            shift[i] = shf / cosphi

        # Now, do it!
        translation = np.zeros(3)
        for i in axes:
            nowlen = lengths[i]
            if vacuum is not None:
                self.cell[i] = cell[i] * (1 + longer[i] / nowlen)
            translation += shift[i] * cell[i] / nowlen

            # We calculated translations using the completed cell,
            # so directions without cell vectors will have been centered
            # along a "fake" vector of length 1.
            # Therefore, we adjust by -0.5:
            if not any(self.cell[i]):
                translation[i] -= 0.5

        # Optionally, translate to center about a point in space.
        if about is not None:
            for vector in self.cell:
                translation -= vector / 2.0
            translation += about

        self.positions += translation

    def get_centroid(self, scaled=False):
        r"""
        Get the centroid of the positions.

        Parameters
        ----------
        scaled : bool
            If scaled=True the centroid in scaled coordinates is returned.

        Notes
        -----
        For a set of :math:`k` positions
        :math:`\textbf{x}_1, \textbf{x}_2, ..., \textbf{x}_k`
        the centroid is given by

        .. math::

            \frac{\textbf{x}_1 + \textbf{x}_2 + ... + \textbf{x}_k}{k}.
        """
        if scaled:
            return self.cell.scaled_positions(self.positions.mean(0))
        return self.positions.mean(0)

    def set_centroid(self, centroid, scaled=False):
        r"""
        Set the centroid of the positions.

        Parameters
        ----------
        centroid : float | np.ndarray
            The new centroid. Can be a float or a xyz vector
        scaled : bool
            If scaled=True the centroid is expected in scaled coordinates.

        Notes
        -----
        For a set of :math:`k` positions
        :math:`\textbf{x}_1, \textbf{x}_2, ..., \textbf{x}_k`
        the centroid is given by

        .. math::

            \frac{\textbf{x}_1 + \textbf{x}_2 + ... + \textbf{x}_k}{k}.
        """
        difference = centroid - self.get_centroid(scaled=scaled)
        if scaled:
            self.set_scaled_positions(self.get_scaled_positions() + difference)
        else:
            self.set_positions(self.get_positions() + difference)

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
            c = np.cos(a)
            s = np.sin(a)
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

    def get_dihedral(self, a0, a1, a2, a3, mic=False):
        """Calculate dihedral angle.

        Calculate dihedral angle (in degrees) between the vectors a0->a1
        and a2->a3.

        Use mic=True to use the Minimum Image Convention and calculate the
        angle across periodic boundaries.
        """
        return self.get_dihedrals([[a0, a1, a2, a3]], mic=mic)[0]

    def get_dihedrals(self, indices, mic=False):
        """Calculate dihedral angles.

        Calculate dihedral angles (in degrees) between the list of vectors
        a0->a1 and a2->a3, where a0, a1, a2 and a3 are in each row of indices.

        Use mic=True to use the Minimum Image Convention and calculate the
        angles across periodic boundaries.
        """
        indices = np.array(indices)
        assert indices.shape[1] == 4

        a0s = self.positions[indices[:, 0]]
        a1s = self.positions[indices[:, 1]]
        a2s = self.positions[indices[:, 2]]
        a3s = self.positions[indices[:, 3]]

        # vectors 0->1, 1->2, 2->3
        v0 = a1s - a0s
        v1 = a2s - a1s
        v2 = a3s - a2s

        cell = None
        pbc = None

        if mic:
            cell = self.cell
            pbc = self.pbc

        return get_dihedrals(v0, v1, v2, cell=cell, pbc=pbc)

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

    def set_dihedral(self, a1, a2, a3, a4, angle, mask=None, indices=None):
        """Set the dihedral angle (degrees) between vectors a1->a2 and
        a3->a4 by changing the qbit indexed by a4.

        If mask is not None, all the qbits described in mask
        (read: the entire subgroup) are moved. Alternatively to the mask,
        the indices of the qbits to be rotated can be supplied. If both
        *mask* and *indices* are given, *indices* overwrites *mask*.

        **Important**: If *mask* or *indices* is given and does not contain
        *a4*, *a4* will NOT be moved. In most cases you therefore want
        to include *a4* in *mask*/*indices*.

        Example: the following defines a very crude
        ethane-like molecule and twists one half of it by 30 degrees.

        >>> qbits = Qbits('HHCCHH', [[-1, 1, 0], [-1, -1, 0], [0, 0, 0],
        ...                          [1, 0, 0], [2, 1, 0], [2, -1, 0]])
        >>> qbits.set_dihedral(1, 2, 3, 4, 210, mask=[0, 0, 0, 1, 1, 1])
        """

        angle = _to_rads(angle)

        # if not provided, set mask to the last qbit in the
        # dihedral description
        if mask is None and indices is None:
            mask = np.zeros(len(self))
            mask[a4] = 1
        elif indices is not None:
            mask = [index in indices for index in range(len(self))]

        # compute necessary in dihedral change, from current value
        current = _to_rads(self.get_dihedral(a1, a2, a3, a4))
        diff = angle - current
        axis = self.positions[a3] - self.positions[a2]
        center = self.positions[a3]
        self._masked_rotate(center, axis, diff, mask)

    def rotate_dihedral(self, a1, a2, a3, a4, angle=None, mask=None, indices=None):
        """Rotate dihedral angle.

        Same usage as in :meth:`ase.Qbits.set_dihedral`: Rotate a group by a
        predefined dihedral angle, starting from its current configuration.
        """
        start = self.get_dihedral(a1, a2, a3, a4)
        self.set_dihedral(a1, a2, a3, a4, angle + start, mask, indices)

    def get_angle(self, index_1: int, index_2: int, index_3: int, mic: bool = False):
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
        mic : bool
            Use mic=True to use the Minimum Image Convention and calculate the
            angle across periodic boundaries.

        Notes
        -----
        Let x1, x2, x3 be the vectors describing the positions of the three
        qubits. Then we calcule the angle between x1-x2 and x3-x2.
        """
        return self.get_angles([[index_1, index_2, index_3]], mic=mic)[0]

    def get_angles(self, indices, mic=False):
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

        if mic:
            return get_angles(v12, v32, cell=self.cell, pbc=self.pbc)

        return get_angles(v12, v32, cell=None, pbc=None)

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

        This method adds random displacements to the qbit positions,
        taking a possible constraint into account.  The random numbers are
        drawn from a normal distribution of standard deviation stdev.

        For a parallel calculation, it is important to use the same
        seed on all processors!"""

        if seed is not None and rng is not None:
            raise ValueError("Please do not provide both seed and rng.")

        if rng is None:
            if seed is None:
                seed = 42
            rng = np.random.RandomState(seed)
        positions = self.arrays["positions"]
        self.set_positions(positions + rng.normal(scale=stdev, size=positions.shape))

    def get_distance(self, a0, a1, mic=False, vector=False):
        """Return distance between two qbits.

        Use mic=True to use the Minimum Image Convention.
        vector=True gives the distance vector (from a0 to a1).
        """
        return self.get_distances(a0, [a1], mic=mic, vector=vector)[0]

    def get_distances(self, a, indices, mic=False, vector=False):
        """Return distances of qbit No.i with a list of qbits.

        Use mic=True to use the Minimum Image Convention.
        vector=True gives the distance vector (from a to self[indices]).
        """
        R = self.arrays["positions"]
        p1 = [R[a]]
        p2 = R[indices]

        cell = None
        pbc = None

        if mic:
            cell = self.cell
            pbc = self.pbc

        D, D_len = get_distances(p1, p2, cell=cell, pbc=pbc)

        if vector:
            D.shape = (-1, 3)
            return D
        else:
            D_len.shape = (-1,)
            return D_len

    def get_all_distances(self, mic=False, vector=False):
        """Return distances of all of the qbits with all of the qbits.

        Use mic=True to use the Minimum Image Convention.
        """
        R = self.arrays["positions"]

        cell = None
        pbc = None

        if mic:
            cell = self.cell
            pbc = self.pbc

        D, D_len = get_distances(R, cell=cell, pbc=pbc)

        if vector:
            return D
        else:
            return D_len

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
            D, D_len = find_mic(D, self.cell, self.pbc)
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

    def get_scaled_positions(self, wrap=True):
        """Get positions relative to unit cell.

        If wrap is True, qbits outside the unit cell will be wrapped into
        the cell in those directions with periodic boundary conditions
        so that the scaled coordinates are between zero and one.

        If any cell vectors are zero, the corresponding coordinates
        are evaluated as if the cell were completed using
        ``cell.complete()``.  This means coordinates will be Cartesian
        as long as the non-zero cell vectors span a Cartesian axis or
        plane."""

        fractional = self.cell.scaled_positions(self.positions)

        if wrap:
            for i, periodic in enumerate(self.pbc):
                if periodic:
                    # Yes, we need to do it twice.
                    # See the scaled_positions.py test.
                    fractional[:, i] %= 1.0
                    fractional[:, i] %= 1.0

        return fractional

    def set_scaled_positions(self, scaled):
        """Set positions relative to unit cell."""
        self.positions[:] = self.cell.cartesian_positions(scaled)

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

    # Rajarshi: Removed this for the moment as there is no usage.
    # def get_temperature(self): """Get the temperature in Kelvin."""

    def __eq__(self, other):
        """Check for identity of two qbits objects.

        Identity means: same positions, states, unit cell and
        periodic boundary conditions."""
        if not isinstance(other, Qbits):
            # print("class check")
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

    # @deprecated('Please use qbits.cell.volume')
    # We kind of want to deprecate this, but the ValueError behaviour
    # might be desirable.  Should we do this?
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
        """Set positions directly, bypassing constraints."""
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
