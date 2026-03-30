import numpy as np
import pytest

from qse.signal import Signal, Signals


@pytest.mark.parametrize(
    "values",
    ([2], np.ones(6), [2, -4]),
)
def test_initialization_duration(values):
    """Check that signals are initialized correctly, when only passed values."""
    signal = Signal(values)
    assert np.allclose(signal.values, values)
    assert len(values) == signal.duration


@pytest.mark.parametrize(
    "n_values",
    (1, 2, 3, 10, 100),
)
@pytest.mark.parametrize(
    "duration_per_val",
    (1, 2, 5),
)
def test_initialization(n_values, duration_per_val):
    """Check that signals are initialized correctly."""
    values = np.random.rand(n_values)
    duration = n_values * duration_per_val
    signal = Signal(values, duration)
    assert np.allclose(signal.values, values)
    assert duration == signal.duration
    assert signal.time_per_value == duration_per_val


@pytest.mark.parametrize(
    "values",
    (np.ones(6), np.ones(16)),
)
@pytest.mark.parametrize(
    "duration",
    (2, 8),
)
def test_initialization_fail(values, duration):
    """Check that signals raise errors correctly."""
    with pytest.raises(Exception):
        signal = Signal(values, duration)


@pytest.mark.parametrize(
    "values",
    (np.ones(6), [2, 5, -7, 8]),
)
@pytest.mark.parametrize(
    "duration_per_val",
    (2, 8),
)
@pytest.mark.parametrize(
    "scalar",
    (0.7, 2, -13.67),
)
def test_scalar_operations(values, duration_per_val, scalar):
    """Check scalar operations work with signals."""
    # test addition
    duration = len(values) * duration_per_val
    signal = Signal(values, duration)
    signal = signal + scalar
    assert signal.duration == duration
    assert np.allclose(signal.values - scalar, values)

    # test in-place addition
    signal = Signal(values, duration)
    signal += scalar
    assert signal.duration == duration
    assert np.allclose(signal.values - scalar, values)

    # test multiplication
    signal = Signal(values, duration)
    signal = signal * scalar
    assert signal.duration == duration
    assert np.allclose(signal.values / scalar, values)

    # test in-place multiplication
    signal = Signal(values, duration)
    signal *= scalar
    assert signal.duration == duration
    assert np.allclose(signal.values / scalar, values)
