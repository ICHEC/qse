import numpy as np
import pytest

import qse


def test_pulser():
    """
    Test initializing the Pulser calculator.
    """
    if not qse.calc.pulser.CALCULATOR_AVAILABLE:
        with pytest.raises(Exception, match="Pulser is not installed."):
            qse.calc.Pulser()

    else:
        qse.calc.Pulser(
            amplitude=qse.Signal(np.ones(6)), detuning=qse.Signal(np.ones(6))
        )


def test_myqlm():
    """
    Test initializing the Myqlm calculator.
    """
    if not qse.calc.myqlm.CALCULATOR_AVAILABLE:
        with pytest.raises(Exception, match="myQLM is not installed."):
            qse.calc.Myqlm()

    else:
        qse.calc.Myqlm()
