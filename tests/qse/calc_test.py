import numpy as np
import pytest

import qse


def test_pulser():
    """Test initializing the Pulser calculator."""
    if not qse.calc.pulser.CALCULATOR_AVAILABLE:
        with pytest.raises(Exception, match="Pulser is not installed."):
            qse.calc.Pulser()

    else:
        qse.calc.Pulser()


def test_myqlm():
    """Test initializing the Myqlm calculator."""
    if not qse.calc.myqlm.CALCULATOR_AVAILABLE:
        with pytest.raises(Exception, match="myQLM is not installed."):
            qse.calc.Myqlm()

    else:
        qse.calc.Myqlm()


def _infidelity(state_1, state_2):
    """Calculate the infidelity between two states."""
    return 1.0 - np.abs((np.conj(state_1).T @ state_2).item()) ** 2


@pytest.mark.parametrize("omega", [10.01, 0.453, 5.6])
@pytest.mark.parametrize("delta", [0.12, 5.242, -13.0])
def test_single_qubit(omega, delta):
    """Test the exact calculator on a single qbit."""
    duration = 400

    mag = np.sqrt(delta**2 + omega**2)
    angle = mag * (duration / 1000) * 0.5
    exact_result = (
        np.cos(angle) - 1j * np.sin(angle) * np.array([[delta], [omega]]) / mag
    )

    # Initialise the pulser calculator
    exact_calc = qse.calc.ExactSimulator(
        amplitude=qse.Signal(np.ones(6) * omega, duration),
        detuning=qse.Signal(np.ones(6) * delta, duration),
    )

    # Compute
    exact_calc.calculate()

    assert _infidelity(exact_result, exact_calc.statevector) < 1e-10
