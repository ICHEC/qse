qse.calc.pulser
===============

.. py:module:: qse.calc.pulser

.. autoapi-nested-parse::

   This module Provides QSE calculator interface to Pulser. It is derived from
   the ASE's CP2K calculator, though at the end it may look very different from
   ASE calculator.

   It defines the `Pulser backend <https://pulser.readthedocs.io/en/stable/>`_
   for analog computation.



Classes
-------

.. autoapisummary::

   qse.calc.pulser.Pulser


Module Contents
---------------

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



