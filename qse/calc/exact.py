import numpy as np


class ExactSimulator:
    r"""
    A calculator that evoles the qubit system exactly.
    Only supported for a single qubit.

    Parameters
    ----------
    amplitude: qse.Signal
        The amplitude pulse.
    detuning: qse.Signal
        The detuning pulse.

    Notes
    -----
    For a single qubit the Hamiltonian is

    .. math::
        H = \frac{\Omega}{2} e^{-i\phi}|0\rangle\langle 1| +
        \frac{\Omega}{2} e^{i\phi}|0\rangle\langle 1|
        -\delta|1\rangle\langle 1|

    where :math:`\Omega` is the amplitude, :math:`\phi` the phase and
    :math:`\delta` the detuning.
    Setting the phase to zero and adding a constant :math:`I\delta/2` we get

    .. math::
        H = \frac{\Omega X + \delta Z}{2}

    Then the time evolution unitary is given by

    .. math::
        U(t) = e^{-iH t} =
        \cos(\Delta t/2) I - i \sin(\Delta t/2) \frac{\Omega X + \delta Z}{\Delta}

    Where :math:`\Delta=\sqrt{\delta^2+\Omega^2}`.
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
