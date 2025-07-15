import numpy as np

import qse


def test_pulser_calc():
    """Test the pulser calculator on a chain of qbits."""
    # Define the lattice
    L = 4
    qbits = qse.lattices.chain(4.0, L)

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
    pulser_calc.calculate()
    pulser_calc.get_spins()
    pulser_calc.get_sij()

    assert pulser_calc.spins.shape == (L, 3)
    assert pulser_calc.sij.shape == (L, L)
