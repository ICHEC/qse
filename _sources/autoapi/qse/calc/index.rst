qse.calc
========

.. py:module:: qse.calc

.. autoapi-nested-parse::

   Calculators
   -----------

   Calculators are high level wrappers on different backend
   and let us offload the computational workload.



Submodules
----------

.. toctree::
   :maxdepth: 1

   /autoapi/qse/calc/myqlm/index
   /autoapi/qse/calc/pulser/index


Classes
-------

.. autoapisummary::

   qse.calc.Calculator
   qse.calc.ExactSimulator
   qse.calc.Myqlm
   qse.calc.Pulser
   qse.calc.Qutip


Functions
---------

.. autoapisummary::

   qse.calc.blockade_radius


Package Contents
----------------

.. py:function:: blockade_radius(rabi_frequency, c6=5420158.53)

   Calculate the blockade radius based on the Rabi frequency and C6 coefficient.
   The blockade radius is a measure of the distance at which the Rydberg blockade
   effect occurs, preventing simultaneous excitation of nearby atoms.

   :Parameters: * **rabi_frequency** (*float*) -- The Rabi frequency (in rads/μs) associated with the atomic transition.
                * **c6** (*float, optional*) -- The C6 coefficient (in rads/μs * μm^6) for the van der Waals interaction.
                  Default is 5420158.53.

   :returns: *float* -- The blockade radius (in μm).

   .. rubric:: Notes

   The blockade radius is calculated as:

   .. math::
       r_b = \left( \frac{C_6}{\Omega} \right)^{1/6}

   where :math:`\Omega` is the Rabi frequency and :math:`C_6` is the van der Waals
   coefficient.


.. py:class:: Calculator(qbits, label: str, is_calculator_available: bool, installation_message: str)

   Base-class for all QSE calculators.

   :Parameters: * **label** (*str*) -- Name used for all files.  Not supported by all calculators.
                  May contain a directory, but please use the directory parameter
                  for that instead.
                * **qbits** (*Qbits object*) -- Optional Qbits object to which the calculator will be
                  attached.  When restarting, qbits will get its positions and
                  unit-cell updated from file.


   .. py:method:: calculate()
      :abstractmethod:


      Do the calculation.



   .. py:method:: get_spins()

      Get spin expectation values.
      If the hamiltonian isn't simulated, it triggers simulation first.

      :returns: *np.ndarray* -- Array of Nx3 containing spin expectation values.

      .. seealso:: :obj:`qse.magnetic.get_spins`



   .. py:method:: get_sij()

      Get spin correlation s_ij.
      If the hamiltonian isn't simulated, it triggers simulation first.

      :returns: *np.ndarray* -- Array of NxN shape containing spin correlations.

      .. seealso:: :obj:`qse.magnetic.get_sij`



   .. py:method:: structure_factor_from_sij(L1: int, L2: int, L3: int)

      Get the structure factor.

      :Parameters: * **L1** (*int*) -- Extent of lattice in x direction.
                   * **L2** (*int*) -- Extent of lattice in y direction.
                   * **L3** (*int*) -- Extent of lattice in z direction.

      :returns: *np.ndarray* -- Array containing the structure factor.

      .. seealso:: :obj:`qse.magnetic.structure_factor_from_sij`



.. py:class:: ExactSimulator(amplitude=None, detuning=None)

   A calculator that evoles the qubit system exactly.
   Only supported for a single qubit.

   :Parameters: * **amplitude** (*qse.Signal, qse.Signals*) -- The amplitude pulse.
                * **detuning** (*qse.Signal, qse.Signals*) -- The detuning pulse.

   .. rubric:: Notes

   For a single qubit the Hamiltonian is

   .. math::
       H = \frac{\Omega}{2} e^{-i\phi}|0\rangle\langle 1| +
       \frac{\Omega}{2} e^{i\phi}|0\rangle\langle 1|
       -\delta|1\rangle\langle 1|

   where :math:`\Omega` is the amplitude, :math:`\phi` the phase and
   :math:`\delta` the detuning.
   Setting the phase to zero and adding a constant :math:`I\delta/2` we get

   .. math::
       H = \frac{\Omega X + \delta Z}{2}

   Then the time evolution unitary is given by

   .. math::
       U(t) = e^{-iH t} =
       \cos(\Delta t/2) I - i \sin(\Delta t/2) \frac{\Omega X + \delta Z}{\Delta}

   Where :math:`\Delta=\sqrt{\delta^2+\Omega^2}`.


.. py:class:: Myqlm(qbits=None, amplitude=None, detuning=None, qpu=None, label='myqlm-run', wtimes=True)

   Bases: :py:obj:`qse.calc.calculator.Calculator`


   QSE-Calculator for MyQLM.

   :Parameters: * **qbits** (*qse.Qbits*) -- The qbits object.
                * **amplitude** (*qse.Signal*) -- The amplitude pulse.
                * **detuning** (*qse.Signal*) -- The detuning pulse.
                * **qpu** (*qat.qpus*) -- The Quantum Processing Unit for executing the job.
                * **label** (*str*) -- The label.
                  Defaults to "pulser-run".
                * **wtimes** (*bool*) -- Whether to print the times.
                  Defaults to True.


   .. py:method:: calculate()

      Run the calculation.



.. py:class:: Pulser(qbits=None, amplitude=None, detuning=None, channel='rydberg_global', magnetic_field=None, device=None, emulator=None, label='pulser-run', wtimes=True)

   Bases: :py:obj:`qse.calc.calculator.Calculator`


   QSE-Calculator for pulser.

   :Parameters: * **qbits** (*qse.Qbits*) -- The qbits object.
                * **amplitude** (*qse.Signal, qse.Signals, pulser.waveforms.Waveform*) -- The amplitude pulse.
                * **detuning** (*qse.Signal, qse.Signals, pulser.waveforms.Waveform*) -- The detuning pulse.
                * **channel** (*str*) -- Which channel to use. For example "rydberg_global" for Rydberg or
                  "mw_global" for microwave.
                  Defaults to "rydberg_global".
                * **magnetic_field** (*np.ndarray | list*) -- A magnetic field. Must be a 3-component array or list.
                  Can only be passed when using the Microwave channel ("mw_global").
                * **device** (*pulser.devices.Device*) -- The device.
                  Defaults to pulser.devices.MockDevice.
                * **emulator** -- The emulator.
                  Defaults to QutipEmulator.
                * **label** (*str*) -- The label.
                  Defaults to "pulser-run".
                * **wtimes** (*bool*) -- Whether to print the times.
                  Defaults to True.

   .. rubric:: Examples

   A simple example of using the Pulser calculator:

   .. jupyter-execute::

       import qse
       qbits = qse.lattices.square(4.0, 3, 3)
       duration = 400
       amp = qse.Signal([1.01], duration)
       det = qse.Signal([0.12], duration)
       pcalc = qse.calc.Pulser(
           qbits=qbits,
           amplitude=amp,
           detuning=det)
       pcalc.build_sequence()
       pcalc.calculate()
       pcalc.get_spins()

   .. rubric:: Notes

   Pulser will only use the x-y coordinates of the Qbits object.

   Pulser is an open-source Python software package.
   It provides easy-to-use libraries for designing and
   simulating pulse sequences that act on programmable
   arrays of neutral qbits, a promising platform for
   quantum computation and simulation.
   Online documentation: https://pulser.readthedocs.io
   White paper: Quantum 6, 629 (2022)
   Source code: https://github.com/pasqal-io/Pulser
   License: Apache 2.0
   - see [LICENSE](https://github.com/pasqal-io/Pulser/blob/master/LICENSE)
   for details.


   .. py:method:: build_sequence()

      Build the sequence of operations involving the qubit coordinates,
      amplitude pulse and detuning pulse.



   .. py:method:: calculate(progress=True)

      Run the calculation.



.. py:class:: Qutip(qbits, amplitude, detuning)

   Bases: :py:obj:`qse.calc.calculator.Calculator`


   Base-class for all QSE calculators.

   :Parameters: * **label** (*str*) -- Name used for all files.  Not supported by all calculators.
                  May contain a directory, but please use the directory parameter
                  for that instead.
                * **qbits** (*Qbits object*) -- Optional Qbits object to which the calculator will be
                  attached.  When restarting, qbits will get its positions and
                  unit-cell updated from file.


   .. py:method:: calculate(e_ops=None)

      Do the calculation.



