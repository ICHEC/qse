from math import sqrt

import numpy as np
from ase.geometry import (
    find_mic,
)

__all__ = [
    "FixCartesian",
    "FixBondLength",
    "FixConstraintSingle",
    "FixQbits",
    "FixCom",
    "FixedPlane",
    "FixConstraint",
    "FixedLine",
    "FixBondLengths",
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
