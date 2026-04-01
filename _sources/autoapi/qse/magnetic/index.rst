qse.magnetic
============

.. py:module:: qse.magnetic

.. autoapi-nested-parse::

   Functions for computing magnetic correlation



Functions
---------

.. autoapisummary::

   qse.magnetic.get_basis
   qse.magnetic.get_spins
   qse.magnetic.get_sisj
   qse.magnetic.structure_factor_from_sij
   qse.magnetic.get_number_operator


Module Contents
---------------

.. py:function:: get_basis(nqbits: int, hsize: int = None)

   Returns a boolean array representing basis in qubit product state.
   Originally it is intended for the full qubit space, for which
   hsize = :math:`2^N`, however if we need a subset of the full, then hsize
   can be smaller.

   :Parameters: * **nqbits** (*int*) -- The number of qubits, :math:`N`.
                * **hsize** (*int, optional*) -- The size of the hamiltonian or hilbert space.
                  Defaults to the full Hilbert space, :math:`2^N`.

   :returns: *np.ndarray* -- The basis of shape (hsize, N).


.. py:function:: get_spins(statevector: numpy.ndarray[complex], nqbits: int, ibasis: numpy.ndarray[bool] = None) -> numpy.ndarray[float]

   Get the expectation value of the spin operators :math:`(S_x, S_y, S_z)`.

   :Parameters: * **statevector** (*np.ndarray[complex]*) -- :math:`2^N` sized complex array representing the statevector.
                * **nqbits** (*int*) -- Number of Qubits or Spins, :math:`N`.
                * **ibasis** (*np.ndarray[bool], optional*) -- Boolean array representing product basis for qubits passed for computing.
                  Defaults to the full Hilbert space.

   :returns: *np.ndarray[float]* -- An Nx3 array with expectation values of spin operator.

   .. rubric:: Notes

   With a list of qubits seen as spin 1/2 object, here one calculates the
   expectation value :math:`\langle\psi| (S_x, S_y, S_z) |\psi\rangle`.
   The state is given as follows:
   :math:`|\psi\rangle  = \sum_i \text{statevector}[i] \, \text{ibasis}[i]`.


.. py:function:: get_sisj(statevector: numpy.ndarray[complex], nqbits: int, ibasis: numpy.ndarray[bool] = None) -> numpy.ndarray[float]

   Compute the spin correlation function
   :math:`S_{ij} = \langle\psi| S_i \cdot S_j |\psi\rangle`.

   :Parameters: * **statevector** (*np.ndarray[complex]*) -- :math:`2^N` sized complex array representing the statevector.
                * **nqbits** (*int*) -- Number of Qubits or Spins, :math:`N`.
                * **ibasis** (*np.ndarray[bool], optional*) -- Boolean array representing product basis for qubits passed for computing.
                  Defaults to the full Hilbert space.

   :returns: *np.ndarray[float]* -- An NxN array with computed expectation value of
             :math:`\langle S_i \cdot S_j\rangle`.

   .. rubric:: Notes

   With a list of qubits seen as spin 1/2 objects, the spins operators are
   :math:`S_i = (S_x, S_y, S_z)`. The state :math:`|\psi\rangle` is given as follows:
   :math:`|\psi\rangle  = \sum_i \text{statevector}[i] \, \text{ibasis}[i]`.


.. py:function:: structure_factor_from_sij(L1: int, L2: int, L3: int, qbits: qse.Qbits, s_ij: numpy.ndarray[float]) -> numpy.ndarray[float]

   Computes the structure factor from the spin correlation.
   The structure factor is just fourier transform of the :math:`S_{ij}`
   with:

   .. math::

       S[q] = \frac{1}{N^2} \sum_{ij} S_{ij} \exp{i q \cdot (x_i - x_j)}.

   The (L1, L2, L3) are passed as shape of the lattice, and there is a qubit
   at each of these lattice sites.

   :Parameters: * **L1** (*int*) -- Extent of lattice in x direction.
                * **L2 (int)** -- Extent of lattice in y direction.
                * **L3 (int)** -- Extent of lattice in z direction.
                * **qbits** (*qse.Qbits*) -- The Qbits object representing the lattice.
                * **s_ij** (*np.ndarray[float]*) -- The array with spin correlation.

   :returns: *np.ndarray[float]* -- Returns the structure factor.


.. py:function:: get_number_operator(statevector: numpy.ndarray[complex], nqbits: int, ibasis: numpy.ndarray[bool] = None) -> numpy.ndarray[float]

   Get the expectation value of the number operators.

   :Parameters: * **statevector** (*np.ndarray[complex]*) -- :math:`2^N` sized complex array representing the statevector.
                * **nqbits** (*int*) -- Number of Qubits or Spins, :math:`N`.
                * **ibasis** (*np.ndarray[bool], optional*) -- Boolean array representing product basis for qubits passed for computing.
                  Defaults to the full Hilbert space.

   :returns: *np.ndarray[float]* -- An N array with expectation values of the number operators.

   .. rubric:: Notes

   The number operator for qubit :math:`i` is given by
   :math:`n_i|b_1,...,b_i,...,b_N\rangle=b_i|b_1,...,b_i,...,b_N\rangle`.
   This function returns the vector :math:`\langle\psi|n_i|\psi\rangle`.


