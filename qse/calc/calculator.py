import qse.magnetic as magnetic

# Recognized names of calculators sorted alphabetically:
names = ["myqlm", "pulser", "qiskit"]


special = {
    "myqlm": "MYQLM",
    "pulser": "Pulser",
    "qiskit": "Qiskit",
}


external_calculators = {}


def register_calculator_class(name, cls):
    """Add the class into the database."""
    assert name not in external_calculators
    external_calculators[name] = cls
    names.append(name)
    names.sort()


def get_calculator_class(name):
    """Return calculator class."""
    if name == "pulser":
        from qse.calc.pulser import Pulser as Calculator
    elif name == "myqlm":
        from qse.calc.myqlm import MyQLM as Calculator
    elif name in external_calculators:
        Calculator = external_calculators[name]
    else:
        classname = special.get(name, name.title())
        module = __import__("qse.calculators." + name, {}, None, [classname])
        Calculator = getattr(module, classname)
    return Calculator


class Calculator:
    """
    Base-class for all QSE calculators.

    Parameters
    ----------
    label: str
        Name used for all files.  Not supported by all calculators.
        May contain a directory, but please use the directory parameter
        for that instead.
    qbits: Qbits object
        Optional Qbits object to which the calculator will be
        attached.  When restarting, qbits will get its positions and
        unit-cell updated from file.
    """

    def __init__(
        self,
        qbits,
        label: str,
        is_calculator_available: bool,
        installation_message: str,
    ):
        """ """
        if not is_calculator_available:
            raise Exception(installation_message)

        self._qbits = qbits
        self._label = label

        self.spins = None
        self.sij = None

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, label):
        self._label = label

    @property
    def qbits(self):
        return self._qbits

    @qbits.setter
    def qbits(self, qbits):
        self._qbits = qbits

    def calculate(self):
        """
        Do the calculation.
        """
        raise NotImplementedError

    def get_spins(self):
        """
        Get spin expectation values.
        If the hamiltonian isn't simulated, it triggers simulation first.

        Returns
        -------
        np.ndarray
            Array of Nx3 containing spin expectation values.

        See Also
        --------
        qse.magnetic.get_spins for more details.
        """
        if self.results is None:
            self.calculate()

        return magnetic.get_spins(self.statevector, len(self.qbits))

    def get_sij(self):
        r"""
        Get spin correlation s_ij.
        If the hamiltonian isn't simulated, it triggers simulation first.

        Returns
        -------
        np.ndarray
            Array of NxN shape containing spin correlations.

        See Also
        --------
        qse.magnetic.get_sij for more details.
        """
        if self.results is None:
            self.calculate()

        sij = magnetic.get_sisj(self.statevector, len(self.qbits))
        self.sij = sij  # quick fix. TODO: proper property setup done
        return sij

    def structure_factor_from_sij(self, L1: int, L2: int, L3: int):
        r"""
        Get the structure factor.

        Parameters
        ----------
        L1: int
            Extent of lattice in x direction.
        L2: int
            Extent of lattice in y direction.
        L3: int
            Extent of lattice in z direction.

        Returns
        -------
        np.ndarray
            Array containing the structure factor.

        See Also
        --------
        qse.magnetic.structure_factor_from_sij for more details.
        """
        return magnetic.structure_factor_from_sij(L1, L2, L3, self.qbits, self.sij)
