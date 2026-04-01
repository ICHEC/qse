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



