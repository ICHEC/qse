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
        qbits=qbits,
        amplitude=qse.Signal(np.ones(6) * omega0, duration),
        detuning=qse.Signal(np.ones(6) * delta0, duration),
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

    spins = pulser_calc.get_spins()
    assert spins.shape == (repeats, 3)

    sij = pulser_calc.get_sij()
    assert sij.shape == (repeats, repeats)
