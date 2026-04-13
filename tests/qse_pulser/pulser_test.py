import numpy as np
import pulser
import pytest

import qse


@pytest.mark.parametrize("dim", [1, 2, 3])
@pytest.mark.parametrize("nqbits", [2, 3, 4])
def test_qbits_to_pulser(dim, nqbits):
    """Check the to_pulser method in qbits"""
    positions = np.random.rand(nqbits, dim)
    positions -= positions.mean(0)  # pulser uses centered coords.

    qbits = qse.Qbits(positions)
    reg = qbits.to_pulser()

    assert isinstance(reg, pulser.Register)
    pulser_coords = np.array([i for i in reg.qubits.values()])
    if dim == 1:
        assert np.allclose(
            pulser_coords, np.column_stack([positions, np.zeros(nqbits)])
        )
    elif dim == 2:
        assert np.allclose(pulser_coords, positions)
    else:  # dim == 3
        assert np.allclose(pulser_coords, positions[:, :2])


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
        amplitude=qse.Signal([omega0], duration),
        detuning=qse.Signal([delta0], duration),
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


def test_pulser_calc_qutip():
    """Test the pulser calculator matches the qutip calculator."""
    duration = 400
    omega0 = 10.01
    delta0 = 0.12

    # Define the lattice
    br = qse.calc.blockade_radius(omega0)
    qbits = qse.lattices.square(0.9 * br, 2, 2)

    amp = qse.Signal([omega0], duration)
    det = qse.Signal([delta0], duration)

    qutip_calc = qse.calc.Qutip(amplitude=amp, detuning=det, qbits=qbits)
    result = qutip_calc.calculate()

    # pulser calc
    pulser_calc = qse.calc.Pulser(amplitude=amp, detuning=det, qbits=qbits)
    pulser_calc.build_sequence()
    pulser_calc.calculate()

    overlap = (np.conj(result["state"].flatten()) * pulser_calc.statevector).sum()
    infidelity = 1.0 - np.abs(overlap) ** 2
    assert infidelity < 1e-4


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
    """Test the pulser calculator against the exact calculator for a single qbit."""
    duration = 400
    amplitude = qse.Signal([omega], duration)
    detuning = qse.Signal([delta], duration)

    # pass time in micro seconds.
    exact_calc = qse.calc.ExactSimulator(amplitude=amplitude, detuning=detuning)
    # Compute
    exact_calc.calculate()

    # Define the lattice
    qbits = qse.Qbits(positions=np.zeros((1, 3)))

    # Initialise the pulser calculator
    pulser_calc = qse.calc.Pulser(
        amplitude=amplitude,
        detuning=detuning,
        qbits=qbits,
    )

    # Compute
    pulser_calc.build_sequence()
    pulser_calc.calculate()

    assert _infidelity(exact_calc.statevector, pulser_calc.statevector) < 1e-5


def test_signal_to_pulser():
    """
    Check that a amplitude and detuning signal matches.
    """
    omega_max = 2.0 * 2 * np.pi  # rad/µs
    rabi_frequency = omega_max / 2.0  # rad/µs
    delta_0 = -6 * rabi_frequency  # ns
    delta_f = 2 * rabi_frequency  # ns
    t_rise = 252  # ns
    t_fall = 500  # ns
    t_sweep = (delta_f - delta_0) / (2 * np.pi * 10) * 1000  # ns
    t_sweep = int(t_sweep)

    # up ramp, constant, downramp waveform
    amplitude_afm = pulser.CompositeWaveform(
        pulser.waveforms.RampWaveform(t_rise, 0.0, omega_max),
        pulser.waveforms.ConstantWaveform(t_sweep, omega_max),
        pulser.waveforms.RampWaveform(t_fall, omega_max, 0.0),
    )
    amplitude = qse.Signals()
    amplitude += qse.Signal(omega_max * np.linspace(0, 1, t_rise))
    amplitude += qse.Signal([omega_max], t_sweep)
    amplitude += qse.Signal(omega_max - omega_max * np.linspace(0, 1, t_fall))

    assert amplitude.to_pulser() == amplitude_afm

    # corresponding waveform for detuning
    detuning_afm = pulser.CompositeWaveform(
        pulser.waveforms.ConstantWaveform(t_rise, delta_0),
        pulser.waveforms.RampWaveform(t_sweep, delta_0, delta_f),
        pulser.waveforms.ConstantWaveform(t_fall, delta_f),
    )

    detuning = qse.Signals()
    detuning += qse.Signal([delta_0], t_rise)
    detuning += qse.Signal(delta_0 + (delta_f - delta_0) * np.linspace(0, 1, t_sweep))
    detuning += qse.Signal([delta_f], t_fall)

    assert detuning.to_pulser() == detuning_afm


def test_default_blockade_radius():
    """Check blockade_radius matches the pulser MockDevice for default c6 value."""
    rabi_frequency = 2 * np.pi
    blockade_radius = pulser.devices.MockDevice.rydberg_blockade_radius(rabi_frequency)

    assert np.isclose(blockade_radius, qse.calc.blockade_radius(rabi_frequency))
