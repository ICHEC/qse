qse.calc.messages
=================

.. py:module:: qse.calc.messages

.. autoapi-nested-parse::

   Messages that are thrown when using the calculator.



Exceptions
----------

.. autoapisummary::

   qse.calc.messages.CalculatorError
   qse.calc.messages.CalculatorSetupError
   qse.calc.messages.EnvironmentError
   qse.calc.messages.InputError
   qse.calc.messages.CalculationFailed
   qse.calc.messages.ReadError
   qse.calc.messages.PropertyNotImplementedError
   qse.calc.messages.PropertyNotPresent


Module Contents
---------------

.. py:exception:: CalculatorError

   Bases: :py:obj:`RuntimeError`


   Base class of error types related to QSE calculators.


.. py:exception:: CalculatorSetupError

   Bases: :py:obj:`CalculatorError`


   Calculation cannot be performed with the given parameters.

   Reasons to raise this errors are:
     * The calculator is not properly configured
       (missing executable, environment variables, ...)
     * The given qbits object is not supported
     * Calculator parameters are unsupported

   Typically raised before a calculation.


.. py:exception:: EnvironmentError

   Bases: :py:obj:`CalculatorSetupError`


   Raised if calculator is not properly set up with QSE.
   May be missing an executable or environment variables.


.. py:exception:: InputError

   Bases: :py:obj:`CalculatorSetupError`


   Raised if inputs given to the calculator were incorrect.
   Bad input keywords or values, or missing pseudopotentials.
   This may be raised before or during calculation, depending on
   when the problem is detected.


.. py:exception:: CalculationFailed

   Bases: :py:obj:`CalculatorError`


   Calculation failed unexpectedly.

   Reasons to raise this error are:
     * Calculation did not converge
     * Calculation ran out of memory
     * Segmentation fault or other abnormal termination
     * Arithmetic trouble (singular matrices, NaN, ...)

   Typically raised during calculation.


.. py:exception:: ReadError

   Bases: :py:obj:`CalculatorError`


   Unexpected irrecoverable error while reading calculation results.


.. py:exception:: PropertyNotImplementedError

   Bases: :py:obj:`NotImplementedError`


   Raised if a calculator does not implement the requested property.


.. py:exception:: PropertyNotPresent

   Bases: :py:obj:`CalculatorError`


   Requested property is missing.

   Maybe it was never calculated, or for some reason was not extracted
   with the rest of the results, without being a fatal ReadError.


