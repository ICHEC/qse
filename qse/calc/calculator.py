import copy
import os
import pathlib
import subprocess
import warnings
from typing import Any, Dict, List, Optional, Set

import numpy as np
from ase.outputs import Properties, all_outputs

import qse.magnetic as magnetic
from qse.calc.messages import (
    CalculationFailed,
    CalculatorSetupError,
    PropertyNotImplementedError,
    PropertyNotPresent,
)


def compare_qbits(qbits1, qbits2, tol=1e-15, excluded_properties=None):
    """Check for system changes since last calculation.  Properties in
    ``excluded_properties`` are not checked."""
    if qbits1 is None:
        system_changes = all_changes[:]
    else:
        system_changes = []

        properties_to_check = set(all_changes)
        if excluded_properties:
            properties_to_check -= set(excluded_properties)

        # Check properties that aren't in Qbits.arrays but are attributes of
        # Qbits objects
        for prop in ["cell", "pbc"]:
            if prop in properties_to_check:
                properties_to_check.remove(prop)
                if not equal(getattr(qbits1, prop), getattr(qbits2, prop), atol=tol):
                    system_changes.append(prop)

        arrays1 = set(qbits1.arrays)
        arrays2 = set(qbits2.arrays)

        # Add any properties that are only in qbits1.arrays or only in
        # qbits2.arrays (and aren't excluded).  Note that if, e.g. arrays1 has
        # `initial_charges` which is merely zeros and arrays2 does not have
        # this array, we'll still assume that the system has changed.  However,
        # this should only occur rarely.
        system_changes += properties_to_check & (arrays1 ^ arrays2)

        # Finally, check all of the non-excluded properties shared by the qbits
        # arrays.
        for prop in properties_to_check & arrays1 & arrays2:
            if not equal(qbits1.arrays[prop], qbits2.arrays[prop], atol=tol):
                system_changes.append(prop)

    return system_changes


all_properties = [
    "energy",
    "energies",
]  # Rajarshi: this needs to populate according to usage


all_changes = ["labels", "positions", "states", "cell", "pbc"]

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


def equal(a, b, tol=None, rtol=None, atol=None):
    """ndarray-enabled comparison function."""
    # XXX Known bugs:
    #  * Comparing cell objects (pbc not part of array representation)
    #  * Infinite recursion for cyclic dicts
    #  * Can of worms is open
    if tol is not None:
        msg = "Use `equal(a, b, rtol=..., atol=...)` instead of `tol=...`"
        warnings.warn(msg, DeprecationWarning)
        assert (
            rtol is None and atol is None
        ), "Do not use deprecated `tol` with `atol` and/or `rtol`"
        rtol = tol
        atol = tol

    a_is_dict = isinstance(a, dict)
    b_is_dict = isinstance(b, dict)
    if a_is_dict or b_is_dict:
        # Check that both a and b are dicts
        if not (a_is_dict and b_is_dict):
            return False
        if a.keys() != b.keys():
            return False
        return all(equal(a[key], b[key], rtol=rtol, atol=atol) for key in a)

    if np.shape(a) != np.shape(b):
        return False

    if rtol is None and atol is None:
        return np.array_equal(a, b)

    if rtol is None:
        rtol = 0
    if atol is None:
        atol = 0

    return np.allclose(a, b, rtol=rtol, atol=atol)


class Parameters(dict):
    """Dictionary for parameters.

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


class Calculator:
    """
    Base-class for all QSE calculators, adapted from ASE calculators.

    A calculator must raise PropertyNotImplementedError if asked for a
    property that it can't calculate.  So, if calculation of the
    stress tensor has not been implemented, get_stress(qbits) should
    raise PropertyNotImplementedError.  This can be achieved simply by not
    including the string 'stress' in the list implemented_properties
    which is a class member.  These are the names of the standard
    properties: 'energy', 'forces', 'stress', 'dipole', 'charges',
    'magmom' and 'magmoms'.
    """

    implemented_properties: List[str] = []
    "Properties calculator can handle (energy, forces, ...)"

    default_parameters: Dict[str, Any] = {}
    "Default parameters"

    ignored_changes: Set[str] = set()
    "Properties of Qbits which we ignore for the purposes of cache "
    "invalidation with check_state()."

    discard_results_on_any_change = False
    "Whether we purge the results following any change in the set() method.  "
    "Most (file I/O) calculators will probably want this."

    _deprecated = object()

    def __init__(
        self,
        is_calculator_available: bool,
        installation_message: str,
        **kwargs,
    ):
        """
        Basic calculator implementation.
        label: str
            Name used for all files.  Not supported by all calculators.
            May contain a directory, but please use the directory parameter
            for that instead.
        qbits: Qbits object
            Optional Qbits object to which the calculator will be
            attached.  When restarting, qbits will get its positions and
            unit-cell updated from file.
        """
        if not is_calculator_available:
            raise Exception(installation_message)

        # print(kwargs.keys())
        self._qbits = kwargs.get("qbits")  # copy of qbits object from last calculation
        self._label = kwargs.get("label")
        # self._label = label
        # self.qbits = qbits
        # self.label = kwargs.get('label')
        self.results = {}  # calculated properties (energy, forces, ...)
        self.parameters = None  # calculational parameters
        self.prefix = None
        # print(self.label)
        if self.parameters is None:
            # Use default parameters if they were not read from file:
            self.parameters = self.get_default_parameters()

        if self.qbits is not None:
            self.qbits.calc = self

        # self.set(**kwargs)

        if not hasattr(self, "name"):
            self.name = self.__class__.__name__.lower()

        self.sij = None

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, label):
        self._label = label

    def get_default_parameters(self):
        return Parameters(copy.deepcopy(self.default_parameters))

    def todict(self, skip_default=True):
        defaults = self.get_default_parameters()
        dct = {}
        for key, value in self.parameters.items():
            if hasattr(value, "todict"):
                value = value.todict()
            if skip_default:
                default = defaults.get(key, "_no_default_")
                if default != "_no_default_" and equal(value, default):
                    continue
            dct[key] = value
        return dct

    def reset(self):
        """Clear all information from old calculation."""

        self.qbits = None
        self.results = {}

    def read(self, label):
        """
        Read qbits, parameters and calculated properties from output file.

        Read result from self.label file.  Raise ReadError if the file
        is not there.  If the file is corrupted or contains an error
        message from the calculation, a ReadError should also be
        raised.  In case of succes, these attributes must set:

        qbits: Qbits object
            The state of the qbits from last calculation.
        parameters: Parameters object
            The parameter dictionary.
        results: dict
            Calculated properties like energy and forces.

        The FileIOCalculator.read() method will typically read qbits
        and parameters and get the results dict by calling the
        read_results() method.
        """

        self.set_label(label)

    def get_qbits(self):
        if self.qbits is None:
            # raise ValueError('Calculator has no qbits')
            print("qbits is None")
            qbits = self.qbits
        else:
            qbits = self.qbits.copy()
            qbits.calc = self
        return qbits

    @property
    def qbits(self):
        return self._qbits

    @qbits.setter
    def qbits(self, qbits):
        self._qbits = qbits

    @classmethod
    def read_qbits(cls, restart, **kwargs):
        return cls(restart=restart, label=restart, **kwargs).get_qbits()

    def set(self, **kwargs):
        """
        Set parameters like set(key1=value1, key2=value2, ...).

        A dictionary containing the parameters that have been changed
        is returned.

        Subclasses must implement a set() method that will look at the
        chaneged parameters and decide if a call to reset() is needed.
        If the changed parameters are harmless, like a change in
        verbosity, then there is no need to call reset().

        The special keyword 'parameters' can be used to read
        parameters from a file.
        """

        if "parameters" in kwargs:
            filename = kwargs.pop("parameters")
            parameters = Parameters.read(filename)
            parameters.update(kwargs)
            kwargs = parameters

        changed_parameters = {}

        for key, value in kwargs.items():
            oldvalue = self.parameters.get(key)
            if key not in self.parameters or not equal(value, oldvalue):
                changed_parameters[key] = value
                self.parameters[key] = value

        if self.discard_results_on_any_change and changed_parameters:
            self.reset()
        return changed_parameters

    def check_state(self, qbits, tol=1e-15):
        """Check for any system changes since last calculation."""
        return compare_qbits(
            self.qbits, qbits, tol=tol, excluded_properties=set(self.ignored_changes)
        )

    def get_energy(self, qbits=None, force_consistent=False):
        energy = self.get_property("energy", qbits)
        if force_consistent:
            if "free_energy" not in self.results:
                name = self.__class__.__name__
                # XXX but we don't know why the energy is not there.
                # We should raise PropertyNotPresent.  Discuss
                raise PropertyNotImplementedError(
                    'Force consistent/free energy ("free_energy") '
                    "not provided by {0} calculator".format(name)
                )
            return self.results["free_energy"]
        else:
            return energy

    def get_property(self, name, qbits=None, allow_calculation=True):
        if name not in self.implemented_properties:
            raise PropertyNotImplementedError(
                "{} property not implemented".format(name)
            )

        if qbits is None:
            qbits = self.qbits
            system_changes = []
        else:
            system_changes = self.check_state(qbits)
            if system_changes:
                self.reset()
        if name not in self.results:
            if not allow_calculation:
                return None
            self.calculate(qbits, [name], system_changes)

        if name not in self.results:
            # For some reason the calculator was not able to do what we want,
            # and that is OK.
            raise PropertyNotImplementedError(
                "{} not present in this " "calculation".format(name)
            )

        result = self.results[name]
        if isinstance(result, np.ndarray):
            result = result.copy()
        return result

    def calculation_required(self, qbits, properties):
        assert not isinstance(properties, str)
        system_changes = self.check_state(qbits)
        if system_changes:
            return True
        for name in properties:
            if name not in self.results:
                return True
        return False

    def calculate(self, qbits=None, properties=["energy"], system_changes=all_changes):
        """
        Do the calculation.

        Parameters
        ----------
        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these six: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

        Subclasses need to implement this, but can ignore properties
        and system_changes if they want.  Calculated properties should
        be inserted into results dictionary like shown in this dummy
        example::

            self.results = {'energy': 0.0,
                            'forces': np.zeros((len(qbits), 3)),
                            'stress': np.zeros(6),
                            'dipole': np.zeros(3),
                            'charges': np.zeros(len(qbits)),
                            'magmom': 0.0,
                            'magmoms': np.zeros(len(qbits))}

        The subclass implementation should first call this
        implementation to set the qbits attribute and create any missing
        directories.
        """

        if qbits is not None:
            self.qbits = qbits.copy()
        if not os.path.isdir(self._directory):
            os.makedirs(self._directory)

    def calculate_properties(self, qbits, properties):
        """This method is experimental; currently for internal use."""
        for name in properties:
            if name not in all_outputs:
                raise ValueError(f"No such property: {name}")

        # We ignore system changes for now.
        self.calculate(qbits, properties, system_changes=all_changes)

        props = self.export_properties()

        for name in properties:
            if name not in props:
                raise PropertyNotPresent(name)
        return props

    def export_properties(self):
        return Properties(self.results)

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

        sij = magnetic.get_sisj(self.statevector, self.nqbits)
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


class FileIOCalculator(Calculator):
    """Base class for calculators that write/read input/output files."""

    command: Optional[str] = None
    "Command used to start calculation"

    def __init__(
        self,
        restart=None,
        ignore_bad_restart_file=Calculator._deprecated,
        label=None,
        qbits=None,
        command=None,
        **kwargs,
    ):
        """File-IO calculator.

        command: str
            Command used to start calculation.
        """

        Calculator.__init__(
            self, restart, ignore_bad_restart_file, label, qbits, **kwargs
        )

        if command is not None:
            self.command = command
        else:
            name = "ASE_" + self.name.upper() + "_COMMAND"
            self.command = os.environ.get(name, self.command)

    def calculate(self, qbits=None, properties=["energy"], system_changes=all_changes):
        Calculator.calculate(self, qbits, properties, system_changes)
        self.write_input(self.qbits, properties, system_changes)
        if self.command is None:
            raise CalculatorSetupError(
                "Please set ${} environment variable ".format(
                    "ASE_" + self.name.upper() + "_COMMAND"
                )
                + "or supply the command keyword"
            )
        command = self.command
        if "PREFIX" in command:
            command = command.replace("PREFIX", self.prefix)

        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory)
        except OSError as err:
            # Actually this may never happen with shell=True, since
            # probably the shell launches successfully.  But we soon want
            # to allow calling the subprocess directly, and then this
            # distinction (failed to launch vs failed to run) is useful.
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err

        errorcode = proc.wait()

        if errorcode:
            path = os.path.abspath(self.directory)
            msg = (
                'Calculator "{}" failed with command "{}" failed in '
                "{} with error code {}".format(self.name, command, path, errorcode)
            )
            raise CalculationFailed(msg)

        self.read_results()

    def write_input(self, qbits, properties=None, system_changes=None):
        """Write input file(s).

        Call this method first in subclasses so that directories are
        created automatically."""

        absdir = os.path.abspath(self.directory)
        if absdir != os.curdir and not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    def read_results(self):
        """Read energy, forces, ... from output file(s)."""
        pass
