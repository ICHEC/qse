"""Messages that are thrown when using the calculator."""

class CalculatorError(RuntimeError):
    """Base class of error types related to QSE calculators."""


class CalculatorSetupError(CalculatorError):
    """
    Calculation cannot be performed with the given parameters.

    Reasons to raise this errors are:
      * The calculator is not properly configured
        (missing executable, environment variables, ...)
      * The given qbits object is not supported
      * Calculator parameters are unsupported

    Typically raised before a calculation.
    """


class EnvironmentError(CalculatorSetupError):
    """
    Raised if calculator is not properly set up with QSE.
    May be missing an executable or environment variables.
    """


class InputError(CalculatorSetupError):
    """
    Raised if inputs given to the calculator were incorrect.
    Bad input keywords or values, or missing pseudopotentials.
    This may be raised before or during calculation, depending on
    when the problem is detected.
    """


class CalculationFailed(CalculatorError):
    """
    Calculation failed unexpectedly.

    Reasons to raise this error are:
      * Calculation did not converge
      * Calculation ran out of memory
      * Segmentation fault or other abnormal termination
      * Arithmetic trouble (singular matrices, NaN, ...)

    Typically raised during calculation.
    """


class ReadError(CalculatorError):
    """Unexpected irrecoverable error while reading calculation results."""


class PropertyNotImplementedError(NotImplementedError):
    """Raised if a calculator does not implement the requested property."""


class PropertyNotPresent(CalculatorError):
    """
    Requested property is missing.

    Maybe it was never calculated, or for some reason was not extracted
    with the rest of the results, without being a fatal ReadError.
    """
