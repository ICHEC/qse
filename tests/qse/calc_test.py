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


def single_qubit_evolution(time, delta, omega):
    """Exact evolution for a single qubit (with zero phase)."""
    magnitude = np.sqrt(delta**2 + omega**2)
    delta_b = delta / magnitude
    omega_b = omega / magnitude
    angle = magnitude * time * 0.5
    return np.array(
        [
            [np.cos(angle) - 1j * np.sin(angle) * delta_b],
            [-1j * np.sin(angle) * omega_b],
        ]
    )


@pytest.mark.parametrize("omega", [10.01, 0.453, 5.6])
@pytest.mark.parametrize("delta", [0.12, 5.242, -13.0])
def test_single_qubit(omega, delta):
    r"""Test the exact calculator on a single qbit.

    The Driving Hamiltonian is:

    $$
    H_D = \frac{\Omega}{2} e^{-i\phi}|a\rangle\langle b| +
    \frac{\Omega}{2} e^{i\phi}|b\rangle\langle a|
    -\delta|b\rangle\langle b|
    $$

    Setting phase to zero and adding a constant $I\delta/2$ we get
    $$
    H_D = \frac{\Omega}{2} (|a\rangle\langle b| + |b\rangle\langle a|)
    +\frac{\delta}{2}(|a\rangle\langle a|-|b\rangle\langle b|)
    $$
    Let $|a\rangle=|0\rangle$ and $|b\rangle=|1\rangle$, then
    $$H_D = \frac{\Omega X +\delta Z}{2}$$

    Define $\Delta=\sqrt{\delta^2+\Omega^2}$ and
    $$
    \begin{split}
    \bar{\delta} &= \frac{\delta}{\Delta}\\
    \bar{\Omega} &= \frac{\Omega}{\Delta}\\
    \end{split}
    $$
    Then
    $$
    U(t)=e^{-iH_D t} = \exp(-i [\bar{\delta} Z + \bar{\Omega} X]\Delta t/2) =
    \cos(\Delta t/2) I - i \sin(\Delta t/2) [\bar{\delta} Z + \bar{\Omega} X]
    $$
    and
    $$
    U(t)|0\rangle = [\cos(\Delta t/2) - i \sin(\Delta t/2)\bar{\delta}]|0\rangle
    - i \sin(\Delta t/2)\bar{\Omega} |1\rangle
    $$
    """
    duration = 400

    # pass time in micro seconds.
    exact_result = single_qubit_evolution(duration / 1000, delta, omega)

    # Initialise the pulser calculator
    exact_calc = qse.calc.ExactSimulator(
        amplitude=qse.Signal(np.ones(6) * omega, duration),
        detuning=qse.Signal(np.ones(6) * delta, duration),
    )

    # Compute
    exact_calc.calculate()
    exact_calc_result = exact_calc.statevector

    assert _infidelity(exact_result, exact_calc_result) < 1e-5
