"""
Functions for computing magnetic correlation
"""

import numpy as np

import qse


def get_basis(nqbits: int, hsize: int = None):
    """
    Returns a boolean array representing basis in qubit product state.
    Originally it is intended for the full qubit space, for which
    hsize = :math:`2^N`, however if we need a subset of the full, then hsize
    can be smaller.

    Parameters
    ----------
    nqbits: int
        The number of qubits, :math:`N`.
    hsize: int, optional
        The size of the hamiltonian or hilbert space.
        Defaults to the full Hilbert space, :math:`2^N`.

    Returns
    -------
    np.ndarray
        The basis of shape (hsize, N).
    """
    if hsize is None:
        hsize = 2**nqbits
    ibasis = np.empty((hsize, nqbits), dtype=bool)
    for i in range(hsize):
        ibasis[i, :] = np.fromiter(np.binary_repr(i, nqbits), dtype=int).astype(bool)
    return ibasis


def _sxop(b: np.ndarray[bool], i: int):
    """
    The :math:`S_x` operator on a boolian array b (basically bit flip).

    Parameters
    ----------
    b: np.ndarray[bool]
        The array representing the basis.
    i: int
        The index in b where to apply operation, i.e., b[i].

    Returns
    -------
    np.ndarray[bool]
        The array after applying the :math:`S_x` operation.
    """
    s = b.copy()
    s[i] = ~s[i]
    return s


def _syop(b: np.ndarray[bool], i: int):
    """
    The :math:`S_y` opertor on a boolian array b (this gives a flipped bit
    at ith index and a sign).

    Parameters
    ----------
    b: np.ndarray[bool]
        The array representing the basis.
    i: int
        The index in b where to apply operation, i.e., b[i].

    Returns
    -------
    tuple
        Of the form (s, c): where s is the basis after operation and c is sign.
    """
    s = b.copy()
    c = (-1) ** s[i] * 1j
    s[i] = ~s[i]
    return (s, c)


def _get_index(arr: np.ndarray[bool], val):
    """
    Get index where value matches in the array.
    Results in either single integer, or an array containing indices.

    Parameters
    ----------
    arr: np.ndarray
        The array.
    val: Any
        The value to match.

    Returns
    -------
    Any
        Integer, or an array of indices.
    """
    return np.where((arr == val).all(axis=1))[0][0]


def get_spins(
    statevector: np.ndarray[complex], nqbits: int, ibasis: np.ndarray[bool] = None
) -> np.ndarray[float]:
    r"""
    Get the expectation value of the spin operators :math:`(S_x, S_y, S_z)`.

    Parameters
    ----------
    statevector: np.ndarray[complex]
        :math:`2^N` sized complex array representing the statevector.
    nqbits: int
        Number of Qubits or Spins, :math:`N`.
    ibasis: np.ndarray[bool], optional
        Boolean array representing product basis for qubits passed for computing.
        Defaults to the full Hilbert space.

    Returns
    -------
    np.ndarray[float]
        An Nx3 array with expectation values of spin operator.

    Notes
    -----
    With a list of qubits seen as spin 1/2 object, here one calculates the
    expectation value :math:`\langle\psi| (S_x, S_y, S_z) |\psi\rangle`.
    The state is given as follows:
    :math:`|\psi\rangle  = \sum_i \text{statevector}[i] \, \text{ibasis}[i]`.
    """
    if ibasis is None:
        ibasis = get_basis(nqbits)

    szi = np.zeros(nqbits, dtype=float)
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        prob = (c_alpha * c_alpha.conj()).real
        zi = 1 - 2 * b
        szi += prob * zi

    sxi = np.zeros(nqbits, dtype=complex)
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        states = [_sxop(b, i) for i in range(nqbits)]
        indices = [_get_index(ibasis, s) for s in states]
        ci = statevector[indices]
        sxi += c_alpha * ci.conj()

    syi = np.zeros(nqbits, dtype=complex)
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        out = [_syop(b, i) for i in range(nqbits)]
        states = [i[0] for i in out]
        cc = np.array([i[1] for i in out])
        indices = [_get_index(ibasis, s) for s in states]
        ci = statevector[indices]
        syi += c_alpha * ci.conj() * cc

    spins = np.array([sxi, syi, szi], dtype=complex).T.real
    return spins


def get_sisj(
    statevector: np.ndarray[complex],
    nqbits: int,
    ibasis: np.ndarray[bool] = None,
) -> np.ndarray[float]:
    r"""
    Compute the spin correlation function
    :math:`S_{ij} = \langle\psi| S_i \cdot S_j |\psi\rangle`.

    Parameters
    ----------
    statevector: np.ndarray[complex]
        :math:`2^N` sized complex array representing the statevector.
    nqbits: int
        Number of Qubits or Spins, :math:`N`.
    ibasis: np.ndarray[bool], optional
        Boolean array representing product basis for qubits passed for computing.
        Defaults to the full Hilbert space.

    Returns
    -------
    np.ndarray[float]
        An NxN array with computed expectation value of
        :math:`\langle S_i \cdot S_j\rangle`.

    Notes
    -----
    With a list of qubits seen as spin 1/2 objects, the spins operators are
    :math:`S_i = (S_x, S_y, S_z)`. The state :math:`|\psi\rangle` is given as follows:
    :math:`|\psi\rangle  = \sum_i \text{statevector}[i] \, \text{ibasis}[i]`.
    """
    if ibasis is None:
        ibasis = get_basis(nqbits)

    s_ij = np.zeros((nqbits, nqbits), dtype=float)

    # z-z part of the correlation
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l]
        prob = (c_alpha * c_alpha.conj()).real
        zi = 1 - 2 * b
        zizj = np.outer(zi, zi)
        s_ij += 0.75 * prob * zizj

    # x-x and y-y part of the correlation
    for l, b in enumerate(ibasis):
        c_alpha = statevector[l].conj()  # this is $c_a^*$
        states_ij = np.array(
            [_sxop(_sxop(b, i), j) for i in range(nqbits) for j in range(nqbits)]
        )
        indices = np.array([_get_index(ibasis, s) for s in states_ij]).reshape(
            nqbits, nqbits
        )
        cij = statevector[indices]  # this is $c_a'$
        zi = b.copy()
        zip = ~zi
        zij1 = np.outer(zi, zip)  # ai * (1-aj)
        zij2 = np.outer(zip, zi)  # (1-ai) * aj
        zij = zij1 + zij2  # sum of sigma +- and -+ term.
        tmp = ((c_alpha * cij) * zij).real
        s_ij += 0.5 * tmp
    #
    # normalise it to 1 for structure factor.
    sij *= 4 / 3
    return s_ij


def structure_factor_from_sij(
    L1: int, L2: int, L3: int, qbits: qse.Qbits, s_ij: np.ndarray[float]
) -> np.ndarray[float]:
    r"""
    Computes the structure factor from the spin correlation.
    The structure factor is just fourier transform of the :math:`S_{ij}`
    with:

    .. math::

        S[q] = \frac{1}{N^2} \sum_{ij} S_{ij} \exp{i q \cdot (x_i - x_j)}.

    The (L1, L2, L3) are passed as shape of the lattice, and there is a qubit
    at each of these lattice sites.

    Parameters
    ----------
    L1: int
        Extent of lattice in x direction.
    L2 (int):
        Extent of lattice in y direction.
    L3 (int):
        Extent of lattice in z direction.
    qbits: qse.Qbits
        The Qbits object representing the lattice.
    s_ij: np.ndarray[float]
        The array with spin correlation.

    Returns
    -------
    np.ndarray[float]
        Returns the structure factor.
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


def get_number_operator(
    statevector: np.ndarray[complex], nqbits: int, ibasis: np.ndarray[bool] = None
) -> np.ndarray[float]:
    r"""
    Get the expectation value of the number operators.

    Parameters
    ----------
    statevector: np.ndarray[complex]
        :math:`2^N` sized complex array representing the statevector.
    nqbits: int
        Number of Qubits or Spins, :math:`N`.
    ibasis: np.ndarray[bool], optional
        Boolean array representing product basis for qubits passed for computing.
        Defaults to the full Hilbert space.

    Returns
    -------
    np.ndarray[float]
        An N array with expectation values of the number operators.

    Notes
    -----
    The number operator for qubit :math:`i` is given by
    :math:`n_i|b_1,...,b_i,...,b_N\rangle=b_i|b_1,...,b_i,...,b_N\rangle`.
    This function returns the vector :math:`\langle\psi|n_i|\psi\rangle`.
    """
    if ibasis is None:
        ibasis = get_basis(nqbits)

    probs = (np.conj(statevector) * statevector).real
    return (ibasis * probs[:, None]).sum(0)
