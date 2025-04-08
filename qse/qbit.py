import numpy as np

# functions adapted from ASE's Qbit/Qbits styled objects
names = {
    "label": ("labels", "R"),
    "state": ("states", np.array([0, 1], dtype=complex)),
    "position": ("positions", np.zeros(3)),
}


def qbitproperty(name, doc):
    """Helper function to easily create Qbit attribute property."""

    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.set(name, value)

    def deleter(self):
        self.delete(name)

    return property(getter, setter, deleter, doc)


def xyzproperty(index):
    """Helper function to easily create Qbit XYZ-property."""

    def getter(self):
        return self.position[index]

    def setter(self, value):
        self.position[index] = value

    return property(getter, setter, doc="XYZ"[index] + "-coordinate")


class Qbit:
    """
    Class for representing a single qbit.

    Parameters
    ----------
    label: str or int
        Can be a str or an int label.
    state: list or tuple or np.ndarray
        Quantum state of the qubit
    position:
        Sequence of 3 floats qubit position.
    qbits:
        ...
    index:
        ...

    Notes
    -----
    Typically one can create an a qubit object
    just by
    q = Qbit()
    """

    __slots__ = ["data", "qbits", "index"]

    def __init__(
        self,
        label="X",
        state=(1, 0),
        position=(0, 0, 0),
        qbits=None,
        index=None,
    ):

        self.data = d = {}

        if qbits is None:
            # This qbit is not part of any Qbits object:
            if isinstance(label, str):
                d["label"] = label
            else:
                d["label"] = "X"

            if isinstance(state, np.ndarray):
                d["state"] = state / np.linalg.norm(state)  # normalise
            else:
                t = np.array(state, complex)
                d["state"] = t / np.linalg.norm(t)
                del t

            d["position"] = np.array(position, float)
        self.index = index
        self.qbits = qbits

    @property
    def scaled_position(self):
        pos = self.position
        spos = self.qbits.cell.scaled_positions(pos[np.newaxis])
        return spos[0]

    @scaled_position.setter
    def scaled_position(self, value):
        pos = self.qbits.cell.cartesian_positions(value)
        self.position = pos

    def __repr__(self):
        # s = "Qbit('%s', %s, %s" % (self.label, list(self.position), list(self.state))
        s = "Qbit(label='%s'" % (self.label)
        for name in ["position", "state"]:
            value = self.get_raw(name)
            if value is not None:
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                s += ", %s=%s" % (name, value)
        if self.qbits is None:
            s += ")"
        else:
            s += ", index=%d)" % self.index
        return s

    def cut_reference_to_qbits(self):
        """Cut reference to qbits object."""
        for name in names:
            self.data[name] = self.get_raw(name)
        self.index = None
        self.qbits = None

    def get_raw(self, name):
        """Get name attribute, return None if not explicitly set."""
        if self.qbits is None:
            return self.data[name]

        plural = names[name][0]
        if plural in self.qbits.arrays:
            return self.qbits.arrays[plural][self.index]
        else:
            return None

    def get(self, name):
        """Get name attribute, return default if not explicitly set."""
        value = self.get_raw(name)
        if value is None:
            value = names[name][1]
        return value

    def set(self, name, value):
        """Set name attribute to value."""
        if self.qbits is None:
            assert name in names
            if name == "state":
                self.data[name] = value / np.linalg.norm(value)
            else:
                self.data[name] = value
        else:
            plural, default = names[name]
            if plural in self.qbits.arrays:
                array = self.qbits.arrays[plural]
                if plural == "states":
                    array[self.index] = value / np.linalg.norm(value)
                else:
                    array[self.index] = value
            else:
                default = np.asarray(default)
                array = np.zeros((len(self.atoms),) + default.shape, default.dtype)
                array[self.index] = value
                self.atoms.new_array(plural, array)

    def delete(self, name):
        """Delete name attribute."""
        assert self.atoms is None
        assert name not in ["label", "position", "state"]
        self.data[name] = None

    state = qbitproperty("state", "Quantum state of qubit as 2-column")
    label = qbitproperty("label", "Integer label asigned to qubit")
    position = qbitproperty("position", "XYZ-coordinates")
    x = xyzproperty(0)
    y = xyzproperty(1)
    z = xyzproperty(2)
