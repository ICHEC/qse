qse.calc.myqlm
==============

.. py:module:: qse.calc.myqlm

.. autoapi-nested-parse::

   Myqlm Calculator
   ----------------

   This module Provides QSE calculator interface to MyQLM backend.

   Ref: https://myqlm.github.io/

   `myqlm` is a software middleware to potentially several quantum
   hardware providers, such as Pasqal devices for analog quantum
   computing.

   It lets one offload an analog quantum workflow to myqlm middleware
   which could be simulated ot real quantum device.

   It has been tested in HPC setting for simulated backend, and
   will be tested when we get access to a real Pasqal device.



Classes
-------

.. autoapisummary::

   qse.calc.myqlm.Myqlm


Module Contents
---------------

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



