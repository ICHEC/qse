import numpy as np
import pytest

import qse


@pytest.mark.parametrize(
    "values",
    ([2], np.ones(6), [2, -4]),
)
def test_initialization_duration(values):
    """Check that signals are initialized correctly, when only passed values."""
    signal = qse.Signal(values)
    assert np.allclose(signal.values, values)
    assert len(values) == signal.duration


@pytest.mark.parametrize(
    "values",
    ([2], np.ones(6)),
)
@pytest.mark.parametrize(
    "duration",
    (2, 8),
)
def test_initialization(values, duration):
    """Check that signals are initialized correctly."""
    signal = qse.Signal(values, duration)
    assert np.allclose(signal.values, values)
    assert duration == signal.duration


@pytest.mark.parametrize(
    "values",
    (np.ones(6), [2, 5, -7, 8]),
)
@pytest.mark.parametrize(
    "duration",
    (2, 8),
)
@pytest.mark.parametrize(
    "scalar",
    (0.7, 2, -13.67),
)
def test_scalar_operations(values, duration, scalar):
    """Check scalar operations work with signals."""
    # test addition
    signal = qse.Signal(values, duration)
    signal = signal + scalar
    assert signal.duration == duration
    assert np.allclose(signal.values - scalar, values)

    # test in-place addition
    signal = qse.Signal(values, duration)
    signal += scalar
    assert signal.duration == duration
    assert np.allclose(signal.values - scalar, values)

    # test multiplication
    signal = qse.Signal(values, duration)
    signal = signal * scalar
    assert signal.duration == duration
    assert np.allclose(signal.values / scalar, values)

    # test in-place multiplication
    signal = qse.Signal(values, duration)
    signal *= scalar
    assert signal.duration == duration
    assert np.allclose(signal.values / scalar, values)


@pytest.mark.parametrize(
    "values_1",
    (np.ones(6), [2, 5, -7, 8]),
)
@pytest.mark.parametrize(
    "duration_1",
    (2, 8),
)
@pytest.mark.parametrize(
    "values_2",
    (np.ones(2), [-1, 5, -7, 8], [1.1]),
)
@pytest.mark.parametrize(
    "duration_2",
    (1, 3),
)
def test_signal_addition(values_1, duration_1, values_2, duration_2):
    """Check that signals can be added together."""
    signal_1 = qse.Signal(values_1, duration_1)
    signal_2 = qse.Signal(values_2, duration_2)
    signal = signal_1 + signal_2

    assert signal.duration == duration_1 + duration_2
    assert np.allclose(signal.values, list(values_1) + list(values_2))
