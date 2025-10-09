import pathlib

import qse.magnetic as magnetic


class Parameters(dict):
    """
    Dictionary for parameters.
    Special feature: If param is a Parameters instance, then param.xc
    is a shorthand for param['xc'].
    """

    def __getattr__(self, key):
        if key not in self:
            return dict.__getattribute__(self, key)
        return self[key]

    def __setattr__(self, key, value):
        self[key] = value

    @classmethod
    def read(cls, filename):
        """Read parameters from file."""
        # We use ast to evaluate literals, avoiding eval()
        # for security reasons.
        import ast

        with open(filename) as fd:
            txt = fd.read().strip()
        assert txt.startswith("dict(")
        assert txt.endswith(")")
        txt = txt[5:-1]

        # The tostring() representation "dict(...)" is not actually
        # a literal, so we manually parse that along with the other
        # formatting that we did manually:
        dct = {}
        for line in txt.splitlines():
            key, val = line.split("=", 1)
            key = key.strip()
            val = val.strip()
            if val[-1] == ",":
                val = val[:-1]
            dct[key] = ast.literal_eval(val)

        parameters = cls(dct)
        return parameters

    def tostring(self):
        keys = sorted(self)
        return (
            "dict("
            + ",\n     ".join("{}={!r}".format(key, self[key]) for key in keys)
            + ")\n"
        )

    def write(self, filename):
        pathlib.Path(filename).write_text(self.tostring())


def get_calculator_class(name):
    """Return calculator class."""
    if name == "pulser":
        from qse.calc.pulser import Pulser as Calculator
    elif name == "myqlm":
        from qse.calc.myqlm import MyQLM as Calculator
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
