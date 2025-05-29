qse.calc.pulser
===============

.. py:module:: qse.calc.pulser

.. autoapi-nested-parse::

   This module Provides QSE calculator interface to Pulser. It is derived from
   the ASE's CP2K calculator, though at the end it may look very different from
   ASE calculator.
   https://pulser.readthedocs.io/en/stable/



Classes
-------

.. autoapisummary::

   qse.calc.pulser.Pulser


Module Contents
---------------

.. py:class:: Pulser(qbits=None, amplitude=None, detuning=None, device=None, emulator=None, label='pulser-run', wtimes=True)

   Bases: :py:obj:`qse.calc.calculator.Calculator`


   QSE-Calculator for pulser.

   Pulser is an open-source Python software package.
   It provides easy-to-use libraries for designing and
   simulating pulse sequences that act on programmable
   arrays of neutral qbits, a promising platform for
   quantum computation and simulation.
   Online documentation: https://pulser.readthedocs.io
   White paper: Quantum 6, 629 (2022)
   Source code: https://github.com/pasqal-io/Pulser
   License: Apache 2.0
   - see [LICENSE](https://github.com/pasqal-io/Pulser/blob/master/LICENSE) for details

   > TODO: if there is any config related to OMP_NUM_THREADS etc, place them here.

   Arguments:
   auto_write: bool
       Flag to enable the auto-write mode. If enabled the
       ``write()`` routine is called after every
       calculation, which mimics the behavior of the
       ``FileIOCalculator``. Default is ``False``.
   basis_set: str
       Name of the basis set to be use.
       The default is ``DZVP-MOLOPT-SR-GTH``.
   basis_set_file: str
       Filename of the basis set file.
       Default is ``BASIS_MOLOPT``.
       Set the environment variable $CP2K_DATA_DIR
       to enabled automatic file discovered.
   charge: float
       The total charge of the system.  Default is ``0``.
   cutoff: float
       The cutoff of the finest grid level.  Default is ``400 * Rydberg``.
   debug: bool
       Flag to enable debug mode. This will print all
       communication between the CP2K-shell and the
       CP2K-calculator. Default is ``False``.
   force_eval_method: str
       The method CP2K uses to evaluate energies and forces.
       The default is ``Quickstep``, which is CP2K's
       module for electronic structure methods like DFT.

   max_scf: int
       Maximum number of SCF iteration to be performed for
       one optimization. Default is ``50``.
   print_level: str
       PRINT_LEVEL of global output.
       Possible options are:
       DEBUG Everything is written out, useful for debugging purposes only
       HIGH Lots of output
       LOW Little output
       MEDIUM Quite some output
       SILENT Almost no output
       Default is 'LOW'


   .. py:attribute:: implemented_properties
      :value: ['energy', 'state', 'fidality']


      Properties calculator can handle (energy, forces, ...)


   .. py:attribute:: default_parameters

      Default parameters


   .. py:method:: set(**kwargs)

      Set parameters like set(key1=value1, key2=value2, ...).



   .. py:method:: write(label)

      Write qbits, parameters and calculated results into restart files.



   .. py:method:: read(label)

      Read qbits, parameters and calculated results from restart files.



   .. py:method:: calculate(progress=True)

      Do the calculation.
      # system_changes=all_changes -> check it's relevance.



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



   .. py:method:: structure_factor_from_sij(L1, L2, L3)

      Get the structure factor

      Args:
          L1 (int): Extent of lattice in x direction
          L2 (int): Extent of lattice in y direction
          L3 (int): Extent of lattice in z direction

      Returns:
          np.ndarray: Array containing the structure factor
      See :py.func: `qse.magnetic.structure_factor_from_sij` for more details.



