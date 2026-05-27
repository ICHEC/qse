"""
Results
-------

This module contains classes that store results of quantum simulation.

"""

from abc import ABC, abstractmethod
from typing import Any, Callable, Dict, Generator

import numpy as np


class BaseResult(ABC):
    """
    Universal class for QSE's result storage.

    The most common, core methods are listed here.
    They should exist regardless of the execution path.
    """

    def __init__(self, metadata: dict = None) -> None:
        self.metadata = metadata or {}

    @abstractmethod
    def get_counts(self, shots: int = 1024) -> Dict[str, int]:
        """
        Returns bitstring frequencies. Calculated
        i. via math for sim,
        ii. via physical measurement for hardware
        """
        pass

    @abstractmethod
    def get_expectation(self, observable: Any) -> float:
        """Returns expectation value.
        (Exact for sim, statistical estimation for hardware)
        """
        pass

    def save(self, filepath: str):
        """Shared utility: Write metadata and data etc."""
        pass


class SimResult(BaseResult):
    """
    Result from a classical emulator are stored here. It intends
    to be pure, vendor-agnostic, and gives access to the quantum state.
    """

    def __init__(
        self, statevector_func: Callable, counts_func: Callable, metadata: dict = None
    ):
        """init function"""
        super().__init__(metadata)  # Hand metadata to the base class.
        self._get_statevector = statevector_func
        self._get_counts = counts_func
        self._cached_state = None

    def get_statevector(self) -> np.ndarray:
        """Returns the exact dense complex array of the quantum state."""
        # Implementation: Extract and format the array from myqlm/pulser
        if self._cached_state is None:
            self._cached_state = self._get_statevector()
        return self._cached_state

    def get_counts(self, shots: int = 1024) -> Dict[str, int]:
        return self._get_counts(shots)

    def states(self) -> Generator[np.ndarray, None, None]:
        """
        Yields the statevector at each time step dt.
        Crucial for analyzing analog Hamiltonian evolution over time.
        """
        pass


class HardwareResult(BaseResult):
    """
    Result from a physical QPU (e.g., Pasqal, QuEra).
    Strictly limited to classical readout data and machine metadata.
    """

    def __init__(self, counts_func: Callable, metadata: Dict = None) -> None:
        super().__init__(metadata)
        self._get_counts = counts_func

    def get_counts(self, shots: int = 1024) -> Dict[str, int]:
        # Implementation: Extract the actual
        # physical measurement counts returned by the QPU
        return self._get_counts(shots)

    # def get_expectation(self, observable: Any) -> float:
    # Implementation: Statistical estimation based on the physical shot counts
    # pass

    def get_statevector(self) -> np.ndarray:
        """
        Physical QPUs cannot expose a statevector.
        This method explicitly blocks unphysical requests.
        """
        raise NotImplementedError(
            "PhysicsError: Cannot extract an exact statevector from a physical QPU. "
            "Use get_counts() to sample the physical system."
        )

    def get_machine_metadata(self) -> Dict[str, Any]:
        """Returns queue time, calibration data, and hardware status."""
        pass
