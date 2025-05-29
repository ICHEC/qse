qse.calc.myqlm
==============

.. py:module:: qse.calc.myqlm

.. autoapi-nested-parse::

   This module Provides QSE calculator interface to MyQLM. It is derived from
   the ASE's calculator, though at the end it may look very different from
   ASE calculator.
   https://myqlm.github.io/



Classes
-------

.. autoapisummary::

   qse.calc.myqlm.Myqlm


Module Contents
---------------

.. py:class:: Myqlm(qbits=None, amplitude=None, detuning=None, duration=None, qpu=None, analog=True, system='rydberg', label='myqlm-run', wtimes=True)

   Bases: :py:obj:`qse.calc.calculator.Calculator`


   QSE-Calculator for MyQLM


   .. py:attribute:: implemented_properties
      :value: ['energy', 'state', 'fidality']


      Properties calculator can handle (energy, forces, ...)


   .. py:attribute:: default_parameters

      Default parameters


   .. py:method:: calculate(qbits=None, properties=..., system_changes=...)

      Do the calculation.

      properties: list of str
          List of what needs to be calculated.  Can be any combination
          of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
          and 'magmoms'.
      system_changes: list of str
          List of what has changed since last calculation.  Can be
          any combination of these six: 'positions', 'numbers', 'cell',
          'pbc', 'initial_charges' and 'initial_magmoms'.

      Subclasses need to implement this, but can ignore properties
      and system_changes if they want.  Calculated properties should
      be inserted into results dictionary like shown in this dummy
      example::

          self.results = {'energy': 0.0,
                          'forces': np.zeros((len(qbits), 3)),
                          'stress': np.zeros(6),
                          'dipole': np.zeros(3),
                          'charges': np.zeros(len(qbits)),
                          'magmom': 0.0,
                          'magmoms': np.zeros(len(qbits))}

      The subclass implementation should first call this
      implementation to set the qbits attribute and create any missing
      directories.



   .. py:method:: get_spins()

      Get spin expectation values
      If the hamiltonian isn't simulated, it triggers simulation first.

      Returns:
          np.ndarray: Array of Nx3 containing spin expectation values.
      See :py.func: `qse.magnetic.get_spins` for more details.



   .. py:method:: get_sij()

      Get spin correlation s_ij
      If the hamiltonian isn't simulated, it triggers simulation first.

      Returns:
          np.ndarray: Array of NxN shape containing spin correlations.
      See :py.func: `qse.magnetic.get_sij` for more details.



   .. py:method:: structure_factor_from_sij(L1: int, L2: int, L3: int)

      Get the structure factor

      Args:
          L1 (int): Extent of lattice in x direction
          L2 (int): Extent of lattice in y direction
          L3 (int): Extent of lattice in z direction

      Returns:
          np.ndarray: Array containing the structure factor
      See :py.func: `qse.magnetic.structure_factor_from_sij` for more details.



