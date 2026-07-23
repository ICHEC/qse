import numpy as np
from scipy.linalg import expm


def evolve(hamiltonian, duration, n_samples, initial_state=None):
    r"""
    Simulate the time evolution of a quantum system under a
    time-independent Hamiltonian.

    Parameters
    ----------
    hamiltonian : ndarray
        The time-independent Hamiltonian matrix of the system. Must be square.
    duration : float
        Total duration of the simulation.
    n_samples : int
        Number of time steps to sample the evolution.
    initial_state : ndarray, optional
        Initial quantum state vector.
        If None, defaults to the first basis state |0⟩ = [1, 0, 0, ..., 0].

    Returns
    -------
    (np.ndarray, np.ndarray)
        A tuple containing the array of times sampled at and a
        matrix of quantum state vectors at each time step,
        including the initial state (each row is a state vector
        of dimension n).

    Notes
    -----
    The Hamiltonian is of form

    .. math::
        H = \sum_{ij} h_{ij} |i\rangle\langle j|

    Where :math:`|i\rangle` is the ith basis state.

    Examples
    --------
    >>> H = np.array([[0, 1], [1, 0]])  # Pauli-X Hamiltonian
    >>> _, states = evolve(H, 1.0, 10)
    >>> states.shape
    (11, 2)
    """
    if hamiltonian.ndim != 2 or (hamiltonian.shape[0] != hamiltonian.shape[1]):
        raise Exception("The Hamiltonian must be square.")

    dt = duration / n_samples
    unitary = expm(-1j * hamiltonian * dt)

    dim = hamiltonian.shape[0]
    states = np.zeros((n_samples + 1, dim), dtype=complex)
    if initial_state is None:
        states[0] = np.array([1] + [0] * (dim - 1))
    else:
        states[0] = initial_state
        if initial_state.shape != (dim,):
            raise Exception(
                f"initial_state must have shape ({dim},), "
                f"got {initial_state.shape}."
            )
    for i in range(n_samples):
        states[i + 1] = unitary @ states[i]

    return times, states
