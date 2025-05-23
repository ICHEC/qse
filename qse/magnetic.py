#!/usr/bin/env python3
r"""
Functions for computing magnetic correlation
"""
import numpy as np

import qse


def get_basis(hsize: int, N: int):
    """returns a boolean array representing basis in qubit product state.
        Originally it is intended for the full qubit space, for which
        hsize = $2^N$, however if we need a subset of the full, then hsize
        can be smaller.

    Args:
        hsize (int): size of the hamiltonian or hilbert space.
        N (int): number of qubits

    Returns:
        np.ndarray: of shape (hsize, N)
    """
    ibasis = np.empty((hsize, N), dtype=bool)
    for i in range(hsize):
        ibasis[i, :] = np.fromiter(np.binary_repr(i, N), dtype=int).astype(bool)
    return ibasis


def sxop(b: np.ndarray[bool], i: int):
    """Sx operator on a boolian array b.
        (Basically bit flip)
    Args:
        b (np.ndarray[bool]): array representing the basis
        i (int): index in b where to apply operation, i.e., b[i]

    Returns:
        np.ndarray[bool]: array after applying the Sx operation
    """
    s = b.copy()
    s[i] = ~s[i]
    return s


def syop(b: np.ndarray[bool], i: int):
    """Sy opertor on a boolian array b
    (This gives a flipped bit at `i`th index and a sign)

    Args:
        b (np.ndarray[bool]): array representing the basis
        i (int): index in b where to apply operation, i.e., b[i]

    Returns:
        tuple (s, c): where s is the basis after operation and c is sign
    """
    s = b.copy()
    s[i] = ~s[i]
    c = (-1) ** s[i] * 1j
    return (s, c)


def get_index(arr: np.ndarray, val):
    """Get index where value matches in the array.
        Results in either single integer, or an array containing indices.
    Args:
        arr (np.ndarray): The array
        val (Any): the value to match

    Returns:
        Any: Integer, or an array of indices.
    """
    return np.where((arr == val).all(axis=1))[0][0]


def get_spins(
    statevector: np.ndarray[complex], ibasis: np.ndarray[bool], N: int
) -> np.ndarray[float]:
    r"""Get the expectation value of the spin operators (Sx, Sy, Sz).
        With a list of qubits seen as spin 1/2 object, here one
        calculates the expectation value <Psi| (Sx, Sy, Sz) |Psi>
        The state is given as follows:
        :math:`|\psi>  = \sum_i statevector[i] ibasis[i]`
    Args:
        statevector (np.ndarray[complex]): :math:`2^N` sized complex array representing the statevector.
        ibasis (np.ndarray[bool]): Boolean array representing product basis
        for qubits passed for computing.
        N (int): Number of Qubits or Spins

    Returns:
        spins -> np.ndarray[float]: An Nx3 array with expectation values of spin operator.
    """
    szi = np.zeros(N, dtype=float)
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        prob = (c_alpha * c_alpha.conj()).real
        zi = 1 - 2 * b
        szi += prob * zi
    #
    # print(szi)
    sxi = np.zeros(N, dtype=complex)
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        states = [sxop(b, i) for i in range(N)]
        indices = [np.where((ibasis == s).all(axis=1))[0][0] for s in states]
        ci = statevector[indices]
        sxi += c_alpha * ci.conj()
    # print(sxi.real)
    syi = np.zeros(N, dtype=complex)
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        out = [syop(b, i) for i in range(N)]
        states = [i[0] for i in out]
        cc = np.array([i[1] for i in out])
        # print(cc, states)
        indices = [np.where((ibasis == s).all(axis=1))[0][0] for s in states]
        ci = statevector[indices]
        syi += c_alpha * ci.conj() * cc
    # print(syi)
    spins = np.array([sxi, syi, szi], dtype=complex).T.real
    return spins


# checking the correction of the code
def get_sisj(
    statevector: np.ndarray[complex], ibasis: np.ndarray[bool], N: int
) -> np.ndarray[float]:
    r"""Compute spin correlation function
        :math: `S_ij = <\psi| S_i \cdot S_j |\psi>.`
        With a list of qubits seen as spin 1/2 object, the spins operators are S_i = (Sx, Sy, Sz).
        The state |Psi> given as follows:
        :math: `|\psi>  = \sum_i statevector[i] ibasis[i]`

    Args:
        statevector (np.ndarray[complex]): :math: `2^N` sized complex array representing the statevector.
        ibasis (np.ndarray[bool]): Boolean array representing product basis
        for qubits passed for computing.
        N (int): Number of Qubits or Spins

    Returns:
        s_ij -> np.ndarray[float]: An NxN array with computed expectation value of <S_i.S_j>
    """
    s_ij = np.zeros((N, N), dtype=float)
    # z-z part of the correlation
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        prob = (c_alpha * c_alpha.conj()).real
        zi = 1 - 2 * b
        zizj = np.outer(zi, zi)
        s_ij += prob * zizj
    # x-x and y-y part of the correlation
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        states_ij = np.array([sxop(sxop(b, i), j) for i in range(N) for j in range(N)])
        indices = np.array(
            [np.where((ibasis == s).all(axis=1))[0][0] for s in states_ij]
        ).reshape(N, N)
        cij = statevector[indices].conj()
        zi = 1 - 2 * b
        zip = -zi
        zij1 = np.outer(zip, zi)
        zij2 = np.outer(zi, zip)
        zij = zij1 + zij2
        tmp = ((c_alpha * cij) * zij).real
        np.fill_diagonal(tmp, 0)
        s_ij += 2 * tmp
    return s_ij


#


def structure_factor_from_sij(
    L1: int, L2: int, L3: int, qbits: qse.Qbits, s_ij: np.ndarray[float]
) -> np.ndarray[float]:
    """From spin correlation, compute the `structure factor`
        The structure factor is just fourier transform of the s_ij
        :math: `S[q] = \frac{1}{N^2} \sum_{ij} s_{ij} \exp{i q \cdot (x_i - x_j)}`
        The (L1, L2, L3) are passed as shape of the lattice, and there is a qubit
        at each of these lattice sites.
    Args:
        L1 (int): Extent of lattice in x direction
        L2 (int): Extent of lattice in y direction
        L3 (int): Extent of lattice in z direction
        qbits (qse.Qbits): The Qbits object representing the lattice.
        s_ij (np.ndarray[float]): The array with spin correlation.

    Returns:
        np.ndarray[float]: Returns the structure factor.
    """
    assert L1 * L2 * L3 == qbits.nqbits
    normalize = qbits.nqbits**2
    struc_fac = np.empty((L1, L2, L3), dtype=complex)
    Qi = tuple((i1, i2, i3) for i1 in range(L1) for i2 in range(L2) for i3 in range(L3))
    Ri = Qi
    Rvecs = Ri @ qbits.cell
    Qvecs = Qi @ qbits.cell.reciprocal() / (L1, L2, L3)
    for q, qvec in zip(Qi, Qvecs):
        exri = np.exp(1j * Rvecs @ qvec)
        exrj = exri.conj()
        struc_fac[q] = np.einsum("ij,i,j->", s_ij, exri, exrj)
    struc_fac /= normalize
    struc_fac = struc_fac.real.copy()
    return struc_fac


#
