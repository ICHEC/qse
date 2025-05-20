#!/usr/bin/env python3
"""
Functions for computing magnetic correlation
"""
import numpy as np


def get_basis(hsize, N):
    ibasis = np.empty((hsize, N), dtype=bool)
    for i in range(hsize):
        ibasis[i, :] = np.fromiter(np.binary_repr(i, N), dtype=int).astype(bool)
    return ibasis


def sxop(b: np.ndarray[bool], i):
    s = b.copy()
    s[i] = ~s[i]
    return s


def syop(b: np.ndarray[bool], i):
    s = b.copy()
    s[i] = ~s[i]
    c = (-1) ** s[i] * 1j
    return (s, c)


def get_index(arr, val):
    return np.where((arr == val).all(axis=1))[0][0]


def get_spins(statevector, ibasis, N):
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
    si = np.array([sxi, syi, szi], dtype=complex).T.real
    return si


# checking the correction of the code
def get_sisj(statevector, ibasis, N):
    sij = np.zeros((N, N), dtype=float)
    # z-z part of the correlation
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        prob = (c_alpha * c_alpha.conj()).real
        zi = 1 - 2 * b
        zizj = np.outer(zi, zi)
        sij += prob * zizj
    # x-x and y-y part of the correlation
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        states_ij = np.array(
            [sxop(sxop(b, i), j) for i in range(N) for j in range(N)]
        )  # .reshape(N, N, N)
        indices = np.array(
            [np.where((ibasis == s).all(axis=1))[0][0] for s in states_ij]
        ).reshape(N, N)
        cij = statevector[indices].conj()
        zi = 1 - 2 * b
        zip = -zi
        zij1 = np.outer(zip, zi)
        zij2 = np.outer(zi, zip)
        zij = zij1 + zij2
        # print(((c_alpha * cij) * zij).imag)
        tmp = ((c_alpha * cij) * zij).real
        np.fill_diagonal(tmp, 0)  # (c_alpha * c_alpha.conj()).real
        sij += 2 * tmp
    return sij


#


def structure_factor_from_sij(L1, L2, L3, qbits, sij):
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
        struc_fac[q] = np.einsum("ij,i,j->", sij, exri, exrj)
    struc_fac /= normalize
    struc_fac = struc_fac.real.copy()
    return struc_fac


#
