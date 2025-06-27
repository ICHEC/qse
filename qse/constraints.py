from math import sqrt
from warnings import warn

import numpy as np
from ase.geometry import (
    conditional_find_mic,
    find_mic,
    get_angles,
    get_angles_derivatives,
    get_dihedrals,
    get_dihedrals_derivatives,
    get_distances_derivatives,
)
from ase.stress import full_3x3_to_voigt_6_stress, voigt_6_to_full_3x3_stress
from ase.utils.parsemath import eval_expression

__all__ = [
    "FixCartesian",
    "FixBondLength",
    "FixedMode",
    "FixConstraintSingle",
    "FixQbits",
    "FixScaled",
    "FixCom",
    "FixedPlane",
    "FixConstraint",
    "FixedLine",
    "FixBondLengths",
    "FixInternals",
    "Hookean",
    "FixScaledParametricRelations",
    "FixCartesianParametricRelations",
]


def dict2constraint(dct):
    if dct["name"] not in __all__:
        raise ValueError
    return globals()[dct["name"]](**dct["kwargs"])


def slice2enlist(s, n):
    """Convert a slice object into a list of (new, old) tuples."""
    if isinstance(s, slice):
        return enumerate(range(*s.indices(n)))
    return enumerate(s)


def constrained_indices(atoms, only_include=None):
    """Returns a list of indices for the atoms that are constrained
    by a constraint that is applied.  By setting only_include to a
    specific type of constraint you can make it only look for that
    given constraint.
    """
    indices = []
    for constraint in atoms.constraints:
        if only_include is not None:
            if not isinstance(constraint, only_include):
                continue
        indices.extend(np.array(constraint.get_indices()))
    return np.array(np.unique(indices))


class FixConstraint:
    """Base class for classes that fix one or more atoms in some way."""

    def index_shuffle(self, atoms, ind):
        """Change the indices.

        When the ordering of the atoms in the Atoms object changes,
        this method can be called to shuffle the indices of the
        constraints.

        ind -- List or tuple of indices.

        """
        raise NotImplementedError

    def repeat(self, m, n):
        """basic method to multiply by m, needs to know the length
        of the underlying atoms object for the assignment of
        multiplied constraints to work.
        """
        msg = (
            "Repeat is not compatible with your atoms' constraints."
            " Use atoms.set_constraint() before calling repeat to "
            "remove your constraints."
        )
        raise NotImplementedError(msg)

    def adjust_momenta(self, atoms, momenta):
        """Adjusts momenta in identical manner to forces."""
        self.adjust_forces(atoms, momenta)

    def copy(self):
        return dict2constraint(self.todict().copy())


class FixConstraintSingle(FixConstraint):
    """Base class for classes that fix a single atom."""

    def __init__(self, a):
        self.a = a

    def index_shuffle(self, atoms, ind):
        """The atom index must be stored as self.a."""
        newa = None  # Signal error
        if self.a < 0:
            self.a += len(atoms)
        for new, old in slice2enlist(ind, len(atoms)):
            if old == self.a:
                newa = new
                break
        if newa is None:
            raise IndexError("Constraint not part of slice")
        self.a = newa

    def get_indices(self):
        return [self.a]


class FixQbits(FixConstraint):
    """Constraint object for fixing some chosen qbits."""

    def __init__(self, indices=None, mask=None):
        """Constrain chosen qbits.

        Parameters
        ----------
        indices : list of int
           Indices for those qbits that should be constrained.
        mask : list of bool
           One boolean per qbit indicating if the qbit should be
           constrained or not.

        Examples
        --------
        Fix all label 'X' qbits:

        >>> mask = [s == 'X' for s in qbits.labels]
        >>> c = FixQbits(mask=mask)
        >>> qbits.set_constraint(c)

        Fix all qbits with z-coordinate less than 1.0 Angstrom:

        >>> c = FixQbits(mask=qbits.positions[:, 2] < 1.0)
        >>> qbits.set_constraint(c)
        """

        if indices is None and mask is None:
            raise ValueError('Use "indices" or "mask".')
        if indices is not None and mask is not None:
            raise ValueError('Use only one of "indices" and "mask".')

        if mask is not None:
            indices = np.arange(len(mask))[np.asarray(mask, bool)]
        else:
            # Check for duplicates:
            srt = np.sort(indices)
            if (np.diff(srt) == 0).any():
                raise ValueError(
                    "FixQbits: The indices array contained duplicates. "
                    "Perhaps you wanted to specify a mask instead, but "
                    "forgot the mask= keyword."
                )
        self.index = np.asarray(indices, int)

        if self.index.ndim != 1:
            raise ValueError("Wrong argument to FixQbits class!")

    def get_removed_dof(self, qbits):
        return 3 * len(self.index)

    def adjust_positions(self, qbits, new):
        new[self.index] = qbits.positions[self.index]

    # def adjust_forces(self, atoms, forces):
    #    forces[self.index] = 0.0

    def index_shuffle(self, qbits, ind):
        # See docstring of superclass
        index = []
        for new, old in slice2enlist(ind, len(qbits)):
            if old in self.index:
                index.append(new)
        if len(index) == 0:
            raise IndexError("All indices in FixQbits not part of slice")
        self.index = np.asarray(index, int)

    def get_indices(self):
        return self.index

    def __repr__(self):
        return "FixQbits(indices=%s)" % ints2string(self.index)

    def todict(self):
        return {"name": "FixQbits", "kwargs": {"indices": self.index.tolist()}}

    def repeat(self, m, n):
        i0 = 0
        nqbits = 0
        if isinstance(m, int):
            m = (m, m, m)
        index_new = []
        for m2 in range(m[2]):
            for m1 in range(m[1]):
                for m0 in range(m[0]):
                    i1 = i0 + n
                    index_new += [i + nqbits for i in self.index]
                    i0 = i1
                    nqbits += n
        self.index = np.asarray(index_new, int)
        return self

    def delete_qbits(self, indices, qbits):
        """Removes qbit number ind from the index array, if present.

        Required for removing qbits with existing FixQbits constraints.
        """

        i = np.zeros(qbits, int) - 1
        new = np.delete(np.arange(qbits), indices)
        i[new] = np.arange(len(new))
        index = i[self.index]
        self.index = index[index >= 0]
        if len(self.index) == 0:
            return None
        return self


class FixCom(FixConstraint):
    """Constraint class for fixing the center of mass.

    References

    https://pubs.acs.org/doi/abs/10.1021/jp9722824

    """

    def get_removed_dof(self, atoms):
        return 3

    def adjust_positions(self, atoms, new):
        masses = atoms.get_masses()
        old_cm = atoms.get_center_of_mass()
        new_cm = np.dot(masses, new) / masses.sum()
        d = old_cm - new_cm
        new += d

    # def adjust_forces(self, atoms, forces):
    #     m = atoms.get_masses()
    #     mm = np.tile(m, (3, 1)).T
    #     lb = np.sum(mm * forces, axis=0) / sum(m**2)
    #     forces -= mm * lb

    def todict(self):
        return {"name": "FixCom", "kwargs": {}}


def ints2string(x, threshold=None):
    """Convert ndarray of ints to string."""
    if threshold is None or len(x) <= threshold:
        return str(x.tolist())
    return str(x[:threshold].tolist())[:-1] + ", ...]"


class FixBondLengths(FixConstraint):
    maxiter = 500

    def __init__(self, pairs, tolerance=1e-13, bondlengths=None, iterations=None):
        """iterations:
        Ignored"""
        self.pairs = np.asarray(pairs)
        self.tolerance = tolerance
        self.bondlengths = bondlengths

    def get_removed_dof(self, atoms):
        return len(self.pairs)

    def adjust_positions(self, atoms, new):
        old = atoms.positions
        masses = atoms.get_masses()

        if self.bondlengths is None:
            self.bondlengths = self.initialize_bond_lengths(atoms)

        for i in range(self.maxiter):
            converged = True
            for j, ab in enumerate(self.pairs):
                a = ab[0]
                b = ab[1]
                cd = self.bondlengths[j]
                r0 = old[a] - old[b]
                d0, _ = find_mic(r0, atoms.cell, atoms.pbc)
                d1 = new[a] - new[b] - r0 + d0
                m = 1 / (1 / masses[a] + 1 / masses[b])
                x = 0.5 * (cd**2 - np.dot(d1, d1)) / np.dot(d0, d1)
                if abs(x) > self.tolerance:
                    new[a] += x * m / masses[a] * d0
                    new[b] -= x * m / masses[b] * d0
                    converged = False
            if converged:
                break
        else:
            raise RuntimeError("Did not converge")

    # def adjust_momenta(self, atoms, p):
    #     old = atoms.positions
    #     masses = atoms.get_masses()

    #     if self.bondlengths is None:
    #         self.bondlengths = self.initialize_bond_lengths(atoms)

    #     for i in range(self.maxiter):
    #         converged = True
    #         for j, ab in enumerate(self.pairs):
    #             a = ab[0]
    #             b = ab[1]
    #             cd = self.bondlengths[j]
    #             d = old[a] - old[b]
    #             d, _ = find_mic(d, atoms.cell, atoms.pbc)
    #             dv = p[a] / masses[a] - p[b] / masses[b]
    #             m = 1 / (1 / masses[a] + 1 / masses[b])
    #             x = -np.dot(dv, d) / cd**2
    #             if abs(x) > self.tolerance:
    #                 p[a] += x * m * d
    #                 p[b] -= x * m * d
    #                 converged = False
    #         if converged:
    #             break
    #     else:
    #         raise RuntimeError('Did not converge')

    # def adjust_forces(self, atoms, forces):
    #     self.constraint_forces = -forces
    #     self.adjust_momenta(atoms, forces)
    #     self.constraint_forces += forces

    def initialize_bond_lengths(self, atoms):
        bondlengths = np.zeros(len(self.pairs))

        for i, ab in enumerate(self.pairs):
            bondlengths[i] = atoms.get_distance(ab[0], ab[1], mic=True)

        return bondlengths

    def get_indices(self):
        return np.unique(self.pairs.ravel())

    def todict(self):
        return {
            "name": "FixBondLengths",
            "kwargs": {"pairs": self.pairs.tolist(), "tolerance": self.tolerance},
        }

    def index_shuffle(self, atoms, ind):
        """Shuffle the indices of the two atoms in this constraint"""
        map = np.zeros(len(atoms), int)
        map[ind] = 1
        n = map.sum()
        map[:] = -1
        map[ind] = range(n)
        pairs = map[self.pairs]
        self.pairs = pairs[(pairs != -1).all(1)]
        if len(self.pairs) == 0:
            raise IndexError("Constraint not part of slice")


def FixBondLength(a1, a2):
    """Fix distance between atoms with indices a1 and a2."""
    return FixBondLengths([(a1, a2)])


class FixedMode(FixConstraint):
    """Constrain atoms to move along directions orthogonal to
    a given mode only."""

    def __init__(self, mode):
        self.mode = (np.asarray(mode) / np.sqrt((mode**2).sum())).reshape(-1)

    def get_removed_dof(self, atoms):
        return len(atoms)

    def adjust_positions(self, atoms, newpositions):
        newpositions = newpositions.ravel()
        oldpositions = atoms.positions.ravel()
        step = newpositions - oldpositions
        newpositions -= self.mode * np.dot(step, self.mode)

    def adjust_forces(self, atoms, forces):
        forces = forces.ravel()
        forces -= self.mode * np.dot(forces, self.mode)

    def index_shuffle(self, atoms, ind):
        eps = 1e-12
        mode = self.mode.reshape(-1, 3)
        excluded = np.ones(len(mode), dtype=bool)
        excluded[ind] = False
        if (abs(mode[excluded]) > eps).any():
            raise IndexError("All nonzero parts of mode not in slice")
        self.mode = mode[ind].ravel()

    def get_indices(self):
        # This function will never properly work because it works on all
        # atoms and it has no idea how to tell how many atoms it is
        # attached to.  If it is being used, surely the user knows
        # everything is being constrained.
        return []

    def todict(self):
        return {"name": "FixedMode", "kwargs": {"mode": self.mode.tolist()}}

    def __repr__(self):
        return "FixedMode(%s)" % self.mode.tolist()


class FixedPlane(FixConstraintSingle):
    """Constrain an atom index *a* to move in a given plane only.

    The plane is defined by its normal vector *direction*."""

    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def get_removed_dof(self, atoms):
        return 1

    def adjust_positions(self, atoms, newpositions):
        step = newpositions[self.a] - atoms.positions[self.a]
        newpositions[self.a] -= self.dir * np.dot(step, self.dir)

    def adjust_forces(self, atoms, forces):
        forces[self.a] -= self.dir * np.dot(forces[self.a], self.dir)

    def todict(self):
        return {
            "name": "FixedPlane",
            "kwargs": {"a": self.a, "direction": self.dir.tolist()},
        }

    def __repr__(self):
        return "FixedPlane(%d, %s)" % (self.a, self.dir.tolist())


class FixedLine(FixConstraintSingle):
    """Constrain an atom index *a* to move on a given line only.

    The line is defined by its vector *direction*."""

    def __init__(self, a, direction):
        self.a = a
        self.dir = np.asarray(direction) / sqrt(np.dot(direction, direction))

    def get_removed_dof(self, atoms):
        return 2

    def adjust_positions(self, atoms, newpositions):
        step = newpositions[self.a] - atoms.positions[self.a]
        x = np.dot(step, self.dir)
        newpositions[self.a] = atoms.positions[self.a] + x * self.dir

    def adjust_forces(self, atoms, forces):
        forces[self.a] = self.dir * np.dot(forces[self.a], self.dir)

    def __repr__(self):
        return "FixedLine(%d, %s)" % (self.a, self.dir.tolist())

    def todict(self):
        return {
            "name": "FixedLine",
            "kwargs": {"a": self.a, "direction": self.dir.tolist()},
        }


class FixCartesian(FixConstraintSingle):
    "Fix an atom index *a* in the directions of the cartesian coordinates."

    def __init__(self, a, mask=(1, 1, 1)):
        self.a = a
        self.mask = ~np.asarray(mask, bool)

    def get_removed_dof(self, atoms):
        return 3 - self.mask.sum()

    def adjust_positions(self, atoms, new):
        step = new[self.a] - atoms.positions[self.a]
        step *= self.mask
        new[self.a] = atoms.positions[self.a] + step

    def adjust_forces(self, atoms, forces):
        forces[self.a] *= self.mask

    def __repr__(self):
        return "FixCartesian(a={0}, mask={1})".format(self.a, list(~self.mask))

    def todict(self):
        return {
            "name": "FixCartesian",
            "kwargs": {"a": self.a, "mask": ~self.mask.tolist()},
        }


class FixScaled(FixConstraintSingle):
    "Fix an atom index *a* in the directions of the unit vectors."

    def __init__(self, cell, a, mask=(1, 1, 1)):
        self.cell = np.asarray(cell)
        self.a = a
        self.mask = np.array(mask, bool)

    def get_removed_dof(self, atoms):
        return self.mask.sum()

    def adjust_positions(self, atoms, new):
        scaled_old = atoms.cell.scaled_positions(atoms.positions)
        scaled_new = atoms.cell.scaled_positions(new)
        for n in range(3):
            if self.mask[n]:
                scaled_new[self.a, n] = scaled_old[self.a, n]
        new[self.a] = atoms.cell.cartesian_positions(scaled_new)[self.a]

    def adjust_forces(self, atoms, forces):
        # Forces are covarient to the coordinate transformation,
        # use the inverse transformations
        scaled_forces = atoms.cell.cartesian_positions(forces)
        scaled_forces[self.a] *= -(self.mask - 1)
        forces[self.a] = atoms.cell.scaled_positions(scaled_forces)[self.a]

    def todict(self):
        return {
            "name": "FixScaled",
            "kwargs": {
                "a": self.a,
                "cell": self.cell.tolist(),
                "mask": self.mask.tolist(),
            },
        }

    def __repr__(self):
        return "FixScaled(%s, %d, %s)" % (repr(self.cell), self.a, repr(self.mask))


# TODO: Better interface might be to use dictionaries in place of very
# nested lists/tuples
class FixInternals(FixConstraint):
    """Constraint object for fixing multiple internal coordinates.

    Allows fixing bonds, angles, and dihedrals.
    Please provide angular units in degrees using angles_deg and
    dihedrals_deg.
    """

    def __init__(
        self,
        bonds=None,
        angles=None,
        dihedrals=None,
        angles_deg=None,
        dihedrals_deg=None,
        bondcombos=None,
        anglecombos=None,
        dihedralcombos=None,
        mic=False,
        epsilon=1.0e-7,
    ):
        # deprecate public API using radians; degrees is preferred
        warn_msg = "Please specify {} in degrees using the {} argument."
        if angles:
            warn(FutureWarning(warn_msg.format("angles", "angle_deg")))
            angles = np.asarray(angles)
            angles[:, 0] = angles[:, 0] / np.pi * 180
            angles = angles.tolist()
        else:
            angles = angles_deg
        if dihedrals:
            warn(FutureWarning(warn_msg.format("dihedrals", "dihedrals_deg")))
            dihedrals = np.asarray(dihedrals)
            dihedrals[:, 0] = dihedrals[:, 0] / np.pi * 180
            dihedrals = dihedrals.tolist()
        else:
            dihedrals = dihedrals_deg

        self.bonds = bonds or []
        self.angles = angles or []
        self.dihedrals = dihedrals or []
        self.bondcombos = bondcombos or []
        self.anglecombos = anglecombos or []
        self.dihedralcombos = dihedralcombos or []
        self.mic = mic
        self.epsilon = epsilon

        self.n = (
            len(self.bonds)
            + len(self.angles)
            + len(self.dihedrals)
            + len(self.bondcombos)
            + len(self.anglecombos)
            + len(self.dihedralcombos)
        )

        # Initialize these at run-time:
        self.constraints = []
        self.initialized = False

    def get_removed_dof(self, atoms):
        return self.n

    def initialize(self, atoms):
        if self.initialized:
            return
        masses = np.repeat(atoms.get_masses(), 3)
        cell = None
        pbc = None
        if self.mic:
            cell = atoms.cell
            pbc = atoms.pbc
        self.constraints = []
        for data, make_constr in [
            (self.bonds, self.FixBondLengthAlt),
            (self.angles, self.FixAngle),
            (self.dihedrals, self.FixDihedral),
            (self.bondcombos, self.FixBondCombo),
            (self.anglecombos, self.FixAngleCombo),
            (self.dihedralcombos, self.FixDihedralCombo),
        ]:
            for datum in data:
                constr = make_constr(datum[0], datum[1], masses, cell, pbc)
                self.constraints.append(constr)
        self.initialized = True

    def shuffle_definitions(self, shuffle_dic, internal_type):
        dfns = []  # definitions
        for dfn in internal_type:  # e.g. for bond in self.bonds
            append = True
            new_dfn = [dfn[0], list(dfn[1])]
            for old in dfn[1]:
                if old in shuffle_dic:
                    new_dfn[1][dfn[1].index(old)] = shuffle_dic[old]
                else:
                    append = False
                    break
            if append:
                dfns.append(new_dfn)
        return dfns

    def shuffle_combos(self, shuffle_dic, internal_type):
        dfns = []  # definitions
        for dfn in internal_type:  # e.g. for bondcombo in self.bondcombos
            append = True
            all_indices = [idx[0:-1] for idx in dfn[1]]
            new_dfn = [dfn[0], list(dfn[1])]
            for i, indices in enumerate(all_indices):
                for old in indices:
                    if old in shuffle_dic:
                        new_dfn[1][i][indices.index(old)] = shuffle_dic[old]
                    else:
                        append = False
                        break
                if not append:
                    break
            if append:
                dfns.append(new_dfn)
        return dfns

    def index_shuffle(self, atoms, ind):
        # See docstring of superclass
        self.initialize(atoms)
        shuffle_dic = dict(slice2enlist(ind, len(atoms)))
        shuffle_dic = {old: new for new, old in shuffle_dic.items()}
        self.bonds = self.shuffle_definitions(shuffle_dic, self.bonds)
        self.angles = self.shuffle_definitions(shuffle_dic, self.angles)
        self.dihedrals = self.shuffle_definitions(shuffle_dic, self.dihedrals)
        self.bondcombos = self.shuffle_combos(shuffle_dic, self.bondcombos)
        self.anglecombos = self.shuffle_combos(shuffle_dic, self.anglecombos)
        self.dihedralcombos = self.shuffle_combos(shuffle_dic, self.dihedralcombos)
        self.initialized = False
        self.initialize(atoms)
        if len(self.constraints) == 0:
            raise IndexError("Constraint not part of slice")

    def get_indices(self):
        cons = []
        for dfn in self.bonds + self.dihedrals + self.angles:
            cons.extend(dfn[1])
        for dfn in self.bondcombos + self.anglecombos + self.dihedralcombos:
            for partial_dfn in dfn[1]:
                cons.extend(partial_dfn[0:-1])  # last index is the coefficient
        return list(set(cons))

    def todict(self):
        return {
            "name": "FixInternals",
            "kwargs": {
                "bonds": self.bonds,
                "angles": self.angles,
                "dihedrals": self.dihedrals,
                "bondcombos": self.bondcombos,
                "anglecombos": self.anglecombos,
                "dihedralcombos": self.dihedralcombos,
                "mic": self.mic,
                "epsilon": self.epsilon,
            },
        }

    def adjust_positions(self, atoms, new):
        self.initialize(atoms)
        for constraint in self.constraints:
            constraint.prepare_jacobian(atoms.positions)
        for j in range(50):
            maxerr = 0.0
            for constraint in self.constraints:
                constraint.adjust_positions(atoms.positions, new)
                maxerr = max(abs(constraint.sigma), maxerr)
            if maxerr < self.epsilon:
                return
        raise ValueError("Shake did not converge.")

    def adjust_forces(self, atoms, forces):
        """Project out translations and rotations and all other constraints"""
        self.initialize(atoms)
        positions = atoms.positions
        N = len(forces)
        list2_constraints = list(np.zeros((6, N, 3)))
        tx, ty, tz, rx, ry, rz = list2_constraints

        list_constraints = [r.ravel() for r in list2_constraints]

        tx[:, 0] = 1.0
        ty[:, 1] = 1.0
        tz[:, 2] = 1.0
        ff = forces.ravel()

        # Calculate the center of mass
        center = positions.sum(axis=0) / N

        rx[:, 1] = -(positions[:, 2] - center[2])
        rx[:, 2] = positions[:, 1] - center[1]
        ry[:, 0] = positions[:, 2] - center[2]
        ry[:, 2] = -(positions[:, 0] - center[0])
        rz[:, 0] = -(positions[:, 1] - center[1])
        rz[:, 1] = positions[:, 0] - center[0]

        # Normalizing transl., rotat. constraints
        for r in list2_constraints:
            r /= np.linalg.norm(r.ravel())

        # Add all angle, etc. constraint vectors
        for constraint in self.constraints:
            constraint.prepare_jacobian(positions)
            constraint.adjust_forces(positions, forces)
            list_constraints.insert(0, constraint.jacobian)
        # QR DECOMPOSITION - GRAM SCHMIDT

        list_constraints = [r.ravel() for r in list_constraints]
        aa = np.column_stack(list_constraints)
        (aa, bb) = np.linalg.qr(aa)
        # Projection
        hh = []
        for i, constraint in enumerate(self.constraints):
            hh.append(aa[:, i] * np.row_stack(aa[:, i]))

        txx = aa[:, self.n] * np.row_stack(aa[:, self.n])
        tyy = aa[:, self.n + 1] * np.row_stack(aa[:, self.n + 1])
        tzz = aa[:, self.n + 2] * np.row_stack(aa[:, self.n + 2])
        rxx = aa[:, self.n + 3] * np.row_stack(aa[:, self.n + 3])
        ryy = aa[:, self.n + 4] * np.row_stack(aa[:, self.n + 4])
        rzz = aa[:, self.n + 5] * np.row_stack(aa[:, self.n + 5])
        T = txx + tyy + tzz + rxx + ryy + rzz
        for vec in hh:
            T += vec
        ff = np.dot(T, np.row_stack(ff))
        forces[:, :] -= np.dot(T, np.row_stack(ff)).reshape(-1, 3)

    def __repr__(self):
        constraints = repr(self.constraints)
        return "FixInternals(_copy_init=%s, epsilon=%s)" % (
            constraints,
            repr(self.epsilon),
        )

    def __str__(self):
        return "\n".join([repr(c) for c in self.constraints])

    # Classes for internal use in FixInternals
    class FixInternalsBase:
        """Base class for subclasses of FixInternals."""

        def __init__(self, targetvalue, indices, masses, cell, pbc):
            self.targetvalue = targetvalue  # constant target value
            self.indices = [defin[0:-1] for defin in indices]  # indices, defs
            self.coefs = np.asarray([defin[-1] for defin in indices])  # coefs
            self.masses = masses
            self.jacobian = []  # geometric Jacobian matrix, Wilson B-matrix
            self.sigma = 1.0  # difference between current and target value
            self.projected_force = None  # helps optimizers scan along constr.
            self.cell = cell
            self.pbc = pbc

        def finalize_jacobian(self, pos, n_internals, n, derivs):
            """Populate jacobian with derivatives for `n_internals` defined
            internals. n = 2 (bonds), 3 (angles), 4 (dihedrals)."""
            jacobian = np.zeros((n_internals, *pos.shape))
            for i, idx in enumerate(self.indices):
                for j in range(n):
                    jacobian[i, idx[j]] = derivs[i, j]
            jacobian = jacobian.reshape((n_internals, 3 * len(pos)))
            self.jacobian = self.coefs @ jacobian

        def finalize_positions(self, newpos):
            jacobian = self.jacobian / self.masses
            lamda = -self.sigma / np.dot(jacobian, self.jacobian)
            dnewpos = lamda * jacobian
            newpos += dnewpos.reshape(newpos.shape)

        def adjust_forces(self, positions, forces):
            self.projected_force = np.dot(self.jacobian, forces.ravel())
            self.jacobian /= np.linalg.norm(self.jacobian)

    class FixBondCombo(FixInternalsBase):
        """Constraint subobject for fixing linear combination of bond lengths
        within FixInternals.

        sum_i( coef_i * bond_length_i ) = constant
        """

        def prepare_jacobian(self, pos):
            bondvectors = [pos[k] - pos[h] for h, k in self.indices]
            derivs = get_distances_derivatives(
                bondvectors, cell=self.cell, pbc=self.pbc
            )
            self.finalize_jacobian(pos, len(bondvectors), 2, derivs)

        def adjust_positions(self, oldpos, newpos):
            bondvectors = [newpos[k] - newpos[h] for h, k in self.indices]
            (_,), (dists,) = conditional_find_mic(
                [bondvectors], cell=self.cell, pbc=self.pbc
            )
            value = np.dot(self.coefs, dists)
            self.sigma = value - self.targetvalue
            self.finalize_positions(newpos)

        def __repr__(self):
            return "FixBondCombo({}, {}, {})".format(
                repr(self.targetvalue), self.indices, self.coefs
            )

    class FixBondLengthAlt(FixBondCombo):
        """Constraint subobject for fixing bond length within FixInternals.
        Fix distance between atoms with indices a1, a2."""

        def __init__(self, targetvalue, indices, masses, cell, pbc):
            indices = [list(indices) + [1.0]]  # bond definition with coef 1.
            super().__init__(targetvalue, indices, masses, cell=cell, pbc=pbc)

        def __repr__(self):
            return "FixBondLengthAlt({}, {})".format(self.targetvalue, *self.indices)

    class FixAngleCombo(FixInternalsBase):
        """Constraint subobject for fixing linear combination of angles
        within FixInternals.

        sum_i( coef_i * angle_i ) = constant
        """

        def gather_vectors(self, pos):
            v0 = [pos[h] - pos[k] for h, k, l in self.indices]
            v1 = [pos[l] - pos[k] for h, k, l in self.indices]
            return v0, v1

        def prepare_jacobian(self, pos):
            v0, v1 = self.gather_vectors(pos)
            derivs = get_angles_derivatives(v0, v1, cell=self.cell, pbc=self.pbc)
            self.finalize_jacobian(pos, len(v0), 3, derivs)

        def adjust_positions(self, oldpos, newpos):
            v0, v1 = self.gather_vectors(newpos)
            value = get_angles(v0, v1, cell=self.cell, pbc=self.pbc)
            value = np.dot(self.coefs, value)
            self.sigma = value - self.targetvalue
            self.finalize_positions(newpos)

        def __repr__(self):
            return "FixAngleCombo({}, {}, {})".format(
                self.targetvalue, self.indices, self.coefs
            )

    class FixAngle(FixAngleCombo):
        """Constraint object for fixing an angle within
        FixInternals using the SHAKE algorithm.

        SHAKE convergence is potentially problematic for angles very close to
        0 or 180 degrees as there is a singularity in the Cartesian derivative.
        """

        def __init__(self, targetvalue, indices, masses, cell, pbc):
            """Fix atom movement to construct a constant angle."""
            indices = [list(indices) + [1.0]]  # angle definition with coef 1.
            super().__init__(targetvalue, indices, masses, cell=cell, pbc=pbc)

        def __repr__(self):
            return "FixAngle({}, {})".format(self.targetvalue, *self.indices)

    class FixDihedralCombo(FixInternalsBase):
        """Constraint subobject for fixing linear combination of dihedrals
        within FixInternals.

        sum_i( coef_i * dihedral_i ) = constant
        """

        def gather_vectors(self, pos):
            v0 = [pos[k] - pos[h] for h, k, l, m in self.indices]
            v1 = [pos[l] - pos[k] for h, k, l, m in self.indices]
            v2 = [pos[m] - pos[l] for h, k, l, m in self.indices]
            return v0, v1, v2

        def prepare_jacobian(self, pos):
            v0, v1, v2 = self.gather_vectors(pos)
            derivs = get_dihedrals_derivatives(v0, v1, v2, cell=self.cell, pbc=self.pbc)
            self.finalize_jacobian(pos, len(v0), 4, derivs)

        def adjust_positions(self, oldpos, newpos):
            v0, v1, v2 = self.gather_vectors(newpos)
            value = get_dihedrals(v0, v1, v2, cell=self.cell, pbc=self.pbc)
            value = np.dot(self.coefs, value)
            self.sigma = value - self.targetvalue
            self.finalize_positions(newpos)

        def __repr__(self):
            return "FixDihedralCombo({}, {}, {})".format(
                self.targetvalue, self.indices, self.coefs
            )

    class FixDihedral(FixDihedralCombo):
        """Constraint object for fixing a dihedral angle using
        the SHAKE algorithm. This one allows also other constraints.

        SHAKE convergence is potentially problematic for near-undefined
        dihedral angles (i.e. when one of the two angles a012 or a123
        approaches 0 or 180 degrees).
        """

        def __init__(self, targetvalue, indices, masses, cell, pbc):
            indices = [list(indices) + [1.0]]  # dihedral def. with coef 1.
            super().__init__(targetvalue, indices, masses, cell=cell, pbc=pbc)

        def adjust_positions(self, oldpos, newpos):
            v0, v1, v2 = self.gather_vectors(newpos)
            value = get_dihedrals(v0, v1, v2, cell=self.cell, pbc=self.pbc)
            # apply minimum dihedral difference 'convention': (diff <= 180)
            self.sigma = (value - self.targetvalue + 180) % 360 - 180
            self.finalize_positions(newpos)

        def __repr__(self):
            return "FixDihedral({}, {})".format(self.targetvalue, *self.indices)


class FixParametricRelations(FixConstraint):
    def __init__(
        self,
        indices,
        Jacobian,
        const_shift,
        params=None,
        eps=1e-12,
        use_cell=False,
    ):
        """Constrains the degrees of freedom to act in a reduced parameter space defined by the Jacobian

        These constraints are based off the work in: https://arxiv.org/abs/1908.01610

        The constraints linearly maps the full 3N degrees of freedom, where N is number of active
        lattice vectors/atoms onto a reduced subset of M free parameters, where M <= 3*N. The
        Jacobian matrix and constant shift vector map the full set of degrees of freedom onto the
        reduced parameter space.

        Currently the constraint is set up to handle either atomic positions or lattice vectors
        at one time, but not both. To do both simply add a two constraints for each set. This is
        done to keep the mathematics behind the operations separate.

        It would be possible to extend these constraints to allow non-linear transformations
        if functionality to update the Jacobian at each position update was included. This would
        require passing an update function evaluate it every time adjust_positions is callled.
        This is currently NOT supported, and there are no plans to implement it in the future.

        Args:
            indices (list of int): indices of the constrained atoms
                (if not None or empty then cell_indices must be None or Empty)
            Jacobian (np.ndarray(shape=(3*len(indices), len(params)))): The Jacobian describing
                the parameter space transformation
            const_shift (np.ndarray(shape=(3*len(indices)))): A vector describing the constant term
                in the transformation not accounted for in the Jacobian
            params (list of str): parameters used in the parametric representation
                if None a list is generated based on the shape of the Jacobian
            eps (float): a small number to compare the similarity of numbers and set the precision used
                to generate the constraint expressions
            use_cell (bool): if True then act on the cell object
        """
        self.indices = np.array(indices)
        self.Jacobian = np.array(Jacobian)
        self.const_shift = np.array(const_shift)

        assert self.const_shift.shape[0] == 3 * len(self.indices)
        assert self.Jacobian.shape[0] == 3 * len(self.indices)

        self.eps = eps
        self.use_cell = use_cell

        if params is None:
            params = []
            if self.Jacobian.shape[1] > 0:
                int_fmt_str = (
                    "{:0" + str(int(np.ceil(np.log10(self.Jacobian.shape[1])))) + "d}"
                )
                for param_ind in range(self.Jacobian.shape[1]):
                    params.append("param_" + int_fmt_str.format(param_ind))
        else:
            assert len(params) == self.Jacobian.shape[-1]

        self.params = params

        self.Jacobian_inv = (
            np.linalg.inv(self.Jacobian.T @ self.Jacobian) @ self.Jacobian.T
        )

    @classmethod
    def from_expressions(cls, indices, params, expressions, eps=1e-12, use_cell=False):
        """Converts the expressions into a Jacobian Matrix/const_shift vector and constructs a FixParametricRelations constraint

        The expressions must be a list like object of size 3*N and elements must be ordered as:
        [n_0,i; n_0,j; n_0,k; n_1,i; n_1,j; .... ; n_N-1,i; n_N-1,j; n_N-1,k],
        where i, j, and k are the first, second and third component of the atomic position/lattice
        vector. Currently only linear operations are allowed to be included in the expressions so
        only terms like:
            - const * param_0
            - sqrt[const] * param_1
            - const * param_0 +/- const * param_1 +/- ... +/- const * param_M
        where const is any real number and param_0, param_1, ..., param_M are the parameters passed in
        params, are allowed.

        For example, the fractional atomic position constraints for wurtzite are:
        params = ["z1", "z2"]
        expressions = [
            "1.0/3.0", "2.0/3.0", "z1",
            "2.0/3.0", "1.0/3.0", "0.5 + z1",
            "1.0/3.0", "2.0/3.0", "z2",
            "2.0/3.0", "1.0/3.0", "0.5 + z2",
        ]

        For diamond are:
        params = []
        expressions = [
            "0.0", "0.0", "0.0",
            "0.25", "0.25", "0.25",
        ],

        and for stannite are
        params=["x4", "z4"]
        expressions = [
            "0.0", "0.0", "0.0",
            "0.0", "0.5", "0.5",
            "0.75", "0.25", "0.5",
            "0.25", "0.75", "0.5",
            "x4 + z4", "x4 + z4", "2*x4",
            "x4 - z4", "x4 - z4", "-2*x4",
             "0.0", "-1.0 * (x4 + z4)", "x4 - z4",
             "0.0", "x4 - z4", "-1.0 * (x4 + z4)",
        ]

        Args:
            indices (list of int): indices of the constrained atoms
                (if not None or empty then cell_indices must be None or Empty)
            params (list of str): parameters used in the parametric representation
            expressions (list of str): expressions used to convert from the parametric to the real space
                representation
            eps (float): a small number to compare the similarity of numbers and set the precision used
                to generate the constraint expressions
            use_cell (bool): if True then act on the cell object

        Returns:
            cls(
                indices,
                Jacobian generated from expressions,
                const_shift generated from expressions,
                params,
                eps-12,
                use_cell,
            )
        """
        Jacobian = np.zeros((3 * len(indices), len(params)))
        const_shift = np.zeros(3 * len(indices))

        for expr_ind, expression in enumerate(expressions):
            expression = expression.strip()

            # Convert subtraction to addition
            expression = expression.replace("-", "+(-1.0)*")
            if "+" == expression[0]:
                expression = expression[1:]
            elif "(+" == expression[:2]:
                expression = "(" + expression[2:]

            # Explicitly add leading zeros so when replacing param_1 with 0.0 param_11 does not become 0.01
            int_fmt_str = "{:0" + str(int(np.ceil(np.log10(len(params) + 1)))) + "d}"

            param_dct = dict()
            param_map = dict()

            # Construct a standardized param template for A/B filling
            for param_ind, param in enumerate(params):
                param_str = "param_" + int_fmt_str.format(param_ind)
                param_map[param] = param_str
                param_dct[param_str] = 0.0

            # Replace the parameters according to the map
            # Sort by string length (long to short) to prevent cases like x11 becoming f"{param_map["x1"]}1"
            for param in sorted(params, key=lambda s: -1.0 * len(s)):
                expression = expression.replace(param, param_map[param])

            # Partial linearity check
            for express_sec in expression.split("+"):
                in_sec = [param in express_sec for param in param_dct]
                n_params_in_sec = len(np.where(np.array(in_sec))[0])
                if n_params_in_sec > 1:
                    raise ValueError(
                        "The FixParametricRelations expressions must be linear."
                    )

            const_shift[expr_ind] = float(eval_expression(expression, param_dct))

            for param_ind in range(len(params)):
                param_str = "param_" + int_fmt_str.format(param_ind)
                if param_str not in expression:
                    Jacobian[expr_ind, param_ind] = 0.0
                    continue
                param_dct[param_str] = 1.0
                test_1 = float(eval_expression(expression, param_dct))
                test_1 -= const_shift[expr_ind]
                Jacobian[expr_ind, param_ind] = test_1

                param_dct[param_str] = 2.0
                test_2 = float(eval_expression(expression, param_dct))
                test_2 -= const_shift[expr_ind]
                if abs(test_2 / test_1 - 2.0) > eps:
                    raise ValueError(
                        "The FixParametricRelations expressions must be linear."
                    )
                param_dct[param_str] = 0.0

        args = [
            indices,
            Jacobian,
            const_shift,
            params,
            eps,
            use_cell,
        ]
        if cls is FixScaledParametricRelations:
            args = args[:-1]
        return cls(*args)

    @property
    def expressions(self):
        """Generate the expressions represented by the current self.Jacobian and self.const_shift objects"""
        expressions = []
        per = int(round(-1 * np.log10(self.eps)))
        fmt_str = "{:." + str(per + 1) + "g}"
        for index, shift_val in enumerate(self.const_shift):
            exp = ""
            if (
                np.all(np.abs(self.Jacobian[index]) < self.eps)
                or np.abs(shift_val) > self.eps
            ):
                exp += fmt_str.format(shift_val)

            param_exp = ""
            for param_index, jacob_val in enumerate(self.Jacobian[index]):
                abs_jacob_val = np.round(np.abs(jacob_val), per + 1)
                if abs_jacob_val < self.eps:
                    continue

                param = self.params[param_index]
                if param_exp or exp:
                    if jacob_val > -1.0 * self.eps:
                        param_exp += " + "
                    else:
                        param_exp += " - "
                elif (not exp) and (not param_exp) and (jacob_val < -1.0 * self.eps):
                    param_exp += "-"

                if np.abs(abs_jacob_val - 1.0) <= self.eps:
                    param_exp += "{:s}".format(param)
                else:
                    param_exp += (fmt_str + "*{:s}").format(abs_jacob_val, param)

            exp += param_exp

            expressions.append(exp)
        return np.array(expressions).reshape((-1, 3))

    def todict(self):
        """Create a dictionary representation of the constraint"""
        return {
            "name": type(self).__name__,
            "kwargs": {
                "indices": self.indices,
                "params": self.params,
                "Jacobian": self.Jacobian,
                "const_shift": self.const_shift,
                "eps": self.eps,
                "use_cell": self.use_cell,
            },
        }

    def __repr__(self):
        """The str representation of the constraint"""
        if len(self.indices) > 1:
            indices_str = "[{:d}, ..., {:d}]".format(self.indices[0], self.indices[-1])
        else:
            indices_str = "[{:d}]".format(self.indices[0])

        if len(self.params) > 1:
            params_str = "[{:s}, ..., {:s}]".format(self.params[0], self.params[-1])
        elif len(self.params) == 1:
            params_str = "[{:s}]".format(self.params[0])
        else:
            params_str = "[]"

        return "{:s}({:s}, {:s}, ..., {:e})".format(
            type(self).__name__, indices_str, params_str, self.eps
        )


class FixScaledParametricRelations(FixParametricRelations):
    def __init__(
        self,
        indices,
        Jacobian,
        const_shift,
        params=None,
        eps=1e-12,
    ):
        """The fractional coordinate version of FixParametricRelations

        All arguments are the same, but since this is for fractional coordinates use_cell is false
        """
        super(FixScaledParametricRelations, self).__init__(
            indices,
            Jacobian,
            const_shift,
            params,
            eps,
            False,
        )

    def adjust_contravariant(self, cell, vecs, B):
        """Adjust the values of a set of vectors that are contravariant with the unit transformation"""
        scaled = cell.scaled_positions(vecs).flatten()
        scaled = self.Jacobian_inv @ (scaled - B)
        scaled = ((self.Jacobian @ scaled) + B).reshape((-1, 3))

        return cell.cartesian_positions(scaled)

    def adjust_positions(self, atoms, positions):
        """Adjust positions of the atoms to match the constraints"""
        positions[self.indices] = self.adjust_contravariant(
            atoms.cell,
            positions[self.indices],
            self.const_shift,
        )
        positions[self.indices] = self.adjust_B(atoms.cell, positions[self.indices])

    def adjust_B(self, cell, positions):
        """Wraps the positions back to the unit cell and adjust B to keep track of this change"""
        fractional = cell.scaled_positions(positions)
        wrapped_fractional = (fractional % 1.0) % 1.0
        self.const_shift += np.round(wrapped_fractional - fractional).flatten()
        return cell.cartesian_positions(wrapped_fractional)

    def adjust_momenta(self, atoms, momenta):
        """Adjust momenta of the atoms to match the constraints"""
        momenta[self.indices] = self.adjust_contravariant(
            atoms.cell,
            momenta[self.indices],
            np.zeros(self.const_shift.shape),
        )

    def adjust_forces(self, atoms, forces):
        """Adjust forces of the atoms to match the constraints"""
        # Forces are coavarient to the coordinate transformation, use the inverse transformations
        cart2frac_jacob = np.zeros(2 * (3 * len(atoms),))
        for i_atom in range(len(atoms)):
            cart2frac_jacob[
                3 * i_atom : 3 * (i_atom + 1), 3 * i_atom : 3 * (i_atom + 1)
            ] = atoms.cell.T

        jacobian = cart2frac_jacob @ self.Jacobian
        jacobian_inv = np.linalg.inv(jacobian.T @ jacobian) @ jacobian.T

        reduced_forces = jacobian.T @ forces.flatten()
        forces[self.indices] = (jacobian_inv.T @ reduced_forces).reshape(-1, 3)

    def todict(self):
        """Create a dictionary representation of the constraint"""
        dct = super(FixScaledParametricRelations, self).todict()
        del dct["kwargs"]["use_cell"]
        return dct


class FixCartesianParametricRelations(FixParametricRelations):
    def __init__(
        self,
        indices,
        Jacobian,
        const_shift,
        params=None,
        eps=1e-12,
        use_cell=False,
    ):
        """The Cartesian coordinate version of FixParametricRelations"""
        super(FixCartesianParametricRelations, self).__init__(
            indices,
            Jacobian,
            const_shift,
            params,
            eps,
            use_cell,
        )

    def adjust_contravariant(self, vecs, B):
        """Adjust the values of a set of vectors that are contravariant with the unit transformation"""
        vecs = self.Jacobian_inv @ (vecs.flatten() - B)
        vecs = ((self.Jacobian @ vecs) + B).reshape((-1, 3))
        return vecs

    def adjust_positions(self, atoms, positions):
        """Adjust positions of the atoms to match the constraints"""
        if self.use_cell:
            return
        positions[self.indices] = self.adjust_contravariant(
            positions[self.indices],
            self.const_shift,
        )

    def adjust_momenta(self, atoms, momenta):
        """Adjust momenta of the atoms to match the constraints"""
        if self.use_cell:
            return
        momenta[self.indices] = self.adjust_contravariant(
            momenta[self.indices],
            np.zeros(self.const_shift.shape),
        )

    def adjust_forces(self, atoms, forces):
        """Adjust forces of the atoms to match the constraints"""
        if self.use_cell:
            return

        forces_reduced = self.Jacobian.T @ forces[self.indices].flatten()
        forces[self.indices] = (self.Jacobian_inv.T @ forces_reduced).reshape(-1, 3)

    def adjust_cell(self, atoms, cell):
        """Adjust the cell of the atoms to match the constraints"""
        if not self.use_cell:
            return
        cell[self.indices] = self.adjust_contravariant(
            cell[self.indices],
            np.zeros(self.const_shift.shape),
        )

    def adjust_stress(self, atoms, stress):
        """Adjust the stress of the atoms to match the constraints"""
        if not self.use_cell:
            return

        stress_3x3 = voigt_6_to_full_3x3_stress(stress)
        stress_reduced = self.Jacobian.T @ stress_3x3[self.indices].flatten()
        stress_3x3[self.indices] = (self.Jacobian_inv.T @ stress_reduced).reshape(-1, 3)

        stress[:] = full_3x3_to_voigt_6_stress(stress_3x3)


class Hookean(FixConstraint):
    """Applies a Hookean restorative force between a pair of atoms, an atom
    and a point, or an atom and a plane."""

    def __init__(self, a1, a2, k, rt=None):
        """Forces two atoms to stay close together by applying no force if
        they are below a threshold length, rt, and applying a Hookean
        restorative force when the distance between them exceeds rt. Can
        also be used to tether an atom to a fixed point in space or to a
        distance above a plane.

        a1 : int
           Index of atom 1
        a2 : one of three options
           1) index of atom 2
           2) a fixed point in cartesian space to which to tether a1
           3) a plane given as (A, B, C, D) in A x + B y + C z + D = 0.
        k : float
           Hooke's law (spring) constant to apply when distance
           exceeds threshold_length. Units of eV A^-2.
        rt : float
           The threshold length below which there is no force. The
           length is 1) between two atoms, 2) between atom and point.
           This argument is not supplied in case 3. Units of A.

        If a plane is specified, the Hooke's law force is applied if the atom
        is on the normal side of the plane. For instance, the plane with
        (A, B, C, D) = (0, 0, 1, -7) defines a plane in the xy plane with a z
        intercept of +7 and a normal vector pointing in the +z direction.
        If the atom has z > 7, then a downward force would be applied of
        k * (atom.z - 7). The same plane with the normal vector pointing in
        the -z direction would be given by (A, B, C, D) = (0, 0, -1, 7).
        """

        if isinstance(a2, int):
            self._type = "two atoms"
            self.indices = [a1, a2]
        elif len(a2) == 3:
            self._type = "point"
            self.index = a1
            self.origin = np.array(a2)
        elif len(a2) == 4:
            self._type = "plane"
            self.index = a1
            self.plane = a2
        else:
            raise RuntimeError("Unknown type for a2")
        self.threshold = rt
        self.spring = k

    def get_removed_dof(self, atoms):
        return 0

    def todict(self):
        dct = {"name": "Hookean"}
        dct["kwargs"] = {"rt": self.threshold, "k": self.spring}
        if self._type == "two atoms":
            dct["kwargs"]["a1"] = self.indices[0]
            dct["kwargs"]["a2"] = self.indices[1]
        elif self._type == "point":
            dct["kwargs"]["a1"] = self.index
            dct["kwargs"]["a2"] = self.origin
        elif self._type == "plane":
            dct["kwargs"]["a1"] = self.index
            dct["kwargs"]["a2"] = self.plane
        else:
            raise NotImplementedError("Bad type: %s" % self._type)
        return dct

    def adjust_positions(self, atoms, newpositions):
        pass

    def adjust_momenta(self, atoms, momenta):
        pass

    def adjust_forces(self, atoms, forces):
        positions = atoms.positions
        if self._type == "plane":
            A, B, C, D = self.plane
            x, y, z = positions[self.index]
            d = (A * x + B * y + C * z + D) / np.sqrt(A**2 + B**2 + C**2)
            if d < 0:
                return
            magnitude = self.spring * d
            direction = -np.array((A, B, C)) / np.linalg.norm((A, B, C))
            forces[self.index] += direction * magnitude
            return
        if self._type == "two atoms":
            p1, p2 = positions[self.indices]
        elif self._type == "point":
            p1 = positions[self.index]
            p2 = self.origin
        displace, _ = find_mic(p2 - p1, atoms.cell, atoms.pbc)
        bondlength = np.linalg.norm(displace)
        if bondlength > self.threshold:
            magnitude = self.spring * (bondlength - self.threshold)
            direction = displace / np.linalg.norm(displace)
            if self._type == "two atoms":
                forces[self.indices[0]] += direction * magnitude
                forces[self.indices[1]] -= direction * magnitude
            else:
                forces[self.index] += direction * magnitude

    def adjust_potential_energy(self, atoms):
        """Returns the difference to the potential energy due to an active
        constraint. (That is, the quantity returned is to be added to the
        potential energy.)"""
        positions = atoms.positions
        if self._type == "plane":
            A, B, C, D = self.plane
            x, y, z = positions[self.index]
            d = (A * x + B * y + C * z + D) / np.sqrt(A**2 + B**2 + C**2)
            if d > 0:
                return 0.5 * self.spring * d**2
            else:
                return 0.0
        if self._type == "two atoms":
            p1, p2 = positions[self.indices]
        elif self._type == "point":
            p1 = positions[self.index]
            p2 = self.origin
        displace, _ = find_mic(p2 - p1, atoms.cell, atoms.pbc)
        bondlength = np.linalg.norm(displace)
        if bondlength > self.threshold:
            return 0.5 * self.spring * (bondlength - self.threshold) ** 2
        else:
            return 0.0

    def get_indices(self):
        if self._type == "two atoms":
            return self.indices
        elif self._type == "point":
            return self.index
        elif self._type == "plane":
            return self.index

    def index_shuffle(self, atoms, ind):
        # See docstring of superclass
        if self._type == "two atoms":
            newa = [-1, -1]  # Signal error
            for new, old in slice2enlist(ind, len(atoms)):
                for i, a in enumerate(self.indices):
                    if old == a:
                        newa[i] = new
            if newa[0] == -1 or newa[1] == -1:
                raise IndexError("Constraint not part of slice")
            self.indices = newa
        elif (self._type == "point") or (self._type == "plane"):
            newa = -1  # Signal error
            for new, old in slice2enlist(ind, len(atoms)):
                if old == self.index:
                    newa = new
                    break
            if newa == -1:
                raise IndexError("Constraint not part of slice")
            self.index = newa

    def __repr__(self):
        if self._type == "two atoms":
            return "Hookean(%d, %d)" % tuple(self.indices)
        elif self._type == "point":
            return "Hookean(%d) to cartesian" % self.index
        else:
            return "Hookean(%d) to plane" % self.index
