qse.magnetic
============

.. py:module:: qse.magnetic

.. autoapi-nested-parse::

   Functions for computing magnetic correlation



Functions
---------

.. autoapisummary::

   qse.magnetic.get_basis
   qse.magnetic.sxop
   qse.magnetic.syop
   qse.magnetic.get_index
   qse.magnetic.get_spins
   qse.magnetic.get_sisj
   qse.magnetic.structure_factor_from_sij


Module Contents
---------------

.. py:function:: get_basis(hsize: int, N: int)

   returns a boolean array representing basis in qubit product state.
       Originally it is intended for the full qubit space, for which
       hsize = $2^N$, however if we need a subset of the full, then hsize
       can be smaller.

   Args:
       hsize (int): size of the hamiltonian or hilbert space.
       N (int): number of qubits

   Returns:
       np.ndarray: of shape (hsize, N)


.. py:function:: sxop(b: numpy.ndarray[bool], i: int)

   Sx operator on a boolian array b.
       (Basically bit flip)
   Args:
       b (np.ndarray[bool]): array representing the basis
       i (int): index in b where to apply operation, i.e., b[i]

   Returns:
       np.ndarray[bool]: array after applying the Sx operation


.. py:function:: syop(b: numpy.ndarray[bool], i: int)

   Sy opertor on a boolian array b
   (This gives a flipped bit at `i`th index and a sign)

   Args:
       b (np.ndarray[bool]): array representing the basis
       i (int): index in b where to apply operation, i.e., b[i]

   Returns:
       tuple (s, c): where s is the basis after operation and c is sign


.. py:function:: get_index(arr: numpy.ndarray, val)

   Get index where value matches in the array.
       Results in either single integer, or an array containing indices.
   Args:
       arr (np.ndarray): The array
       val (Any): the value to match

   Returns:
       Any: Integer, or an array of indices.


.. py:function:: get_spins(statevector: numpy.ndarray[complex], ibasis: numpy.ndarray[bool], N: int) -> numpy.ndarray[float]

   Get the expectation value of the spin operators (Sx, Sy, Sz).
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


.. py:function:: get_sisj(statevector: numpy.ndarray[complex], ibasis: numpy.ndarray[bool], N: int) -> numpy.ndarray[float]

   Compute spin correlation function
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


.. py:function:: structure_factor_from_sij(L1: int, L2: int, L3: int, qbits: qse.Qbits, s_ij: numpy.ndarray[float]) -> numpy.ndarray[float]

   From spin correlation, compute the `structure factor`.
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


