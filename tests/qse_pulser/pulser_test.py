import numpy as np
import pytest

import qse


@pytest.mark.parametrize("duration", [2, 17])
@pytest.mark.parametrize("signal_type", ["amplitude", "detuning"])
def test_pulser_signals(duration, signal_type):
    """Check initializing and updating the amplitude and detuning signals."""
    signal = qse.Signal(np.arange(duration))

    # check initializing
    pulser_calc = qse.calc.Pulser(**{signal_type: signal})
    assert getattr(pulser_calc, signal_type).duration == duration
    assert np.allclose(getattr(pulser_calc, signal_type).samples, signal.values)

    # check we can update the sigals
    signal_new = qse.Signal(0.7 * np.arange(duration + 2))
    setattr(pulser_calc, signal_type, signal_new)
    assert getattr(pulser_calc, signal_type).duration == signal_new.duration
    assert np.allclose(getattr(pulser_calc, signal_type).samples, signal_new.values)


def test_pulser_calc():
    """Test the pulser calculator on a chain of qbits."""
    # Define the lattice
    repeats = 4
    qbits = qse.lattices.chain(4.0, repeats)

    duration = 400
    omega0 = 10.01
    delta0 = 0.12

    # Initialise the pulser calculator
    pulser_calc = qse.calc.Pulser(
        amplitude=qse.Signal(np.ones(6) * omega0, duration),
        detuning=qse.Signal(np.ones(6) * delta0, duration),
        qbits=qbits,
        label="test_run",
    )
    assert pulser_calc.amplitude.duration == duration
    assert np.allclose(pulser_calc.amplitude.samples[0], omega0)

    assert pulser_calc.detuning.duration == duration
    assert np.allclose(pulser_calc.detuning.samples[0], delta0)

    assert pulser_calc.label == "test_run"

    # Compute
    pulser_calc.build_sequence()
    pulser_calc.calculate()
    pulser_calc.get_spins()
    pulser_calc.get_sij()

    assert pulser_calc.spins.shape == (repeats, 3)
    assert pulser_calc.sij.shape == (repeats, repeats)


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
    r"""Test the pulser calculator on a single qbit.

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

    # Define the lattice
    qbits = qse.Qbits(positions=np.zeros((1, 3)))

    # Initialise the pulser calculator
    pulser_calc = qse.calc.Pulser(
        amplitude=qse.Signal(np.ones(6) * omega, duration),
        detuning=qse.Signal(np.ones(6) * delta, duration),
        qbits=qbits,
    )

    # Compute
    pulser_calc.build_sequence()
    pulser_calc.calculate()
    pulser_result = pulser_calc.statevector

    assert _infidelity(exact_result, pulser_result) < 1e-5
