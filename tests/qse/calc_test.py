import numpy as np
import pytest

import qse


def test_pulser():
    """
    Test initializing the Pulser calculator.
    """
    signal = qse.Signal(np.ones(6))
    if not qse.calc.pulser.CALCULATOR_AVAILABLE:
        with pytest.raises(Exception, match="Pulser is not installed."):
            qse.calc.Pulser(amplitude=signal, detuning=signal, qbits=qse.Qbits())

    else:
        qse.calc.Pulser(amplitude=signal, detuning=signal, qbits=qse.Qbits())


def test_myqlm():
    """
    Test initializing the Myqlm calculator.
    """
    if not qse.calc.myqlm.CALCULATOR_AVAILABLE:
        with pytest.raises(Exception, match="myQLM is not installed."):
            qse.calc.Myqlm()

    else:
        qse.calc.Myqlm()
