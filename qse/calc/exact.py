import numpy as np


class ExactSimulator:
    """
    A calculator that evoles the qubit system exactly.
    Only supported for a single qubit.

    Parameters
    ----------
    amplitude: qse.Signal
        The amplitude pulse.
    detuning: qse.Signal
        The detuning pulse.
    """

    def __init__(self, amplitude=None, detuning=None):
        _check_pulses(amplitude, detuning)
        self.amplitude = amplitude
        self.detuning = detuning

    def calculate(self):
        delta_t = _nano_to_micro(self.amplitude.duration / len(self.amplitude.values))
        state = np.array([[1.0], [0.0]])
        for amp, det in zip(self.amplitude.values, self.detuning.values):
            unitary = _get_unitary(amp, det, delta_t)
            state = unitary @ state
        self.statevector = state


def _get_unitary(omega, delta, time):
    """Exact unitary for constant delta, omega, over a given time."""
    magnitude = np.sqrt(delta**2 + omega**2)
    angle = magnitude * time * 0.5
    op = np.array([[delta, omega], [omega, -delta]]) / magnitude
    return np.cos(angle) * np.eye(2) - 1j * np.sin(angle) * op


def _check_pulses(amplitude, detuning):
    if not amplitude.duration == detuning.duration:
        raise Exception("The amplitude and detuning must have the same duration.")
    if not len(amplitude.values) == len(detuning.values):
        raise Exception(
            "The amplitude and detuning must have the same amount of values."
        )


def _nano_to_micro(t):
    return t / 1000
