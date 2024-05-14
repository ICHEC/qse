
import numpy as np
# import pulser

# functions adapted from ASE's Qbit/Qbits styled objects
names = {'label': ('labels', 'R'),
         'state': ('states', np.array([0, 1], dtype=complex)),
        'position': ('positions', np.zeros(3)),
         'tag': ('tags', 0)}

def qbitproperty(name, doc):
    """Helper function to easily create Qbit attribute property."""

    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.set(name, value)

    def deleter(self):
        self.delete(name)

    return property(getter, setter, deleter, doc)


def abcproperty(index):
    """Helper function to easily create Qbit ABC-property."""

    def getter(self):
        return self.scaled_position[index]

    def setter(self, value):
        # We can't just do self.scaled_position[i] = value
        # because scaled_position is a new buffer, not a view into
        # something we can write back to.
        # This is a clear bug!
        spos = self.scaled_position
        spos[index] = value
        self.scaled_position = spos

    return property(getter, setter, doc='ABC'[index] + '-coordinate')


def xyzproperty(index):
    """Helper function to easily create Qbit XYZ-property."""

    def getter(self):
        return self.position[index]

    def setter(self, value):
        self.position[index] = value

    return property(getter, setter, doc='XYZ'[index] + '-coordinate')


class Qbit:
    """Class for representing a single qbit.

    Parameters:

    label: str or int
        Can be a str or an int label.
    position: sequence of 3 floats qubit position.
    tag: int
        Special purpose tag.
    Typically one can create an a qubit object
    just by
    q = Qbit()
    
    """
    __slots__ = ['data', 'qbits', 'index']

    def __init__(self, label='X', state=(1, 0), position=(0, 0, 0), tag=None, qbits=None, index=None):

        self.data = d = {}

        if qbits is None:
            # This qbit is not part of any Qbits object:
            if isinstance(label, str):
                d['label'] = label
            else:
                d['label'] = 'X'
            if isinstance(state, np.ndarray):
                d['state'] = state / np.linalg.norm(state) # normalise
            else:
                t = np.array(state, complex)
                d['state'] = t / np.linalg.norm(t)
                del t
            #
            d['position'] = np.array(position, float)
            d['tag'] = tag
            d['label'] = str(label)
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
        #s = "Qbit('%s', %s, %s" % (self.label, list(self.position), list(self.state))
        s = "Qbit(label='%s'" % (self.label)
        for name in ['position', 'state', 'tag']:
            value = self.get_raw(name)
            if value is not None:
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                s += ', %s=%s' % (name, value)
        if self.qbits is None:
            s += ')'
        else:
            s += ', index=%d)' % self.index
        return s

    def cut_reference_to_qbits(self):
        """Cut reference to qbits object."""
        for name in names:
            self.data[name] = self.get_raw(name)
        self.index = None
        self.qbits = None

    def get_raw(self, name):
        """Get name attribute, return None if not explicitly set."""
        #if name == 'symbol':
        #    return chemical_symbols[self.get_raw('number')]

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
        #if name == 'symbol':
        #    name = 'number'
        #    value = atomic_numbers[value]

        if self.qbits is None:
            assert name in names
            if name == 'state':
                self.data[name] = value / np.linalg.norm(value)
            else:
                self.data[name] = value
        else:
            plural, default = names[name]
            if plural in self.qbits.arrays:
                array = self.qbits.arrays[plural]
                #if name == 'magmom' and array.ndim == 2:
                #    assert len(value) == 3
                if plural == 'states':
                    array[self.index] = value / np.linalg.norm(value)
                else:
                    array[self.index] = value
            else:
                #if name == 'magmom' and np.asarray(value).ndim == 1:
                #    array = np.zeros((len(self.atoms), 3))
                #elif name == 'mass':
                #    array = self.atoms.get_masses()
                #else:
                default = np.asarray(default)
                array = np.zeros((len(self.atoms),) + default.shape,
                                     default.dtype)
                array[self.index] = value
                self.atoms.new_array(plural, array)

    def delete(self, name):
        """Delete name attribute."""
        assert self.atoms is None
        assert name not in ['label', 'tag', 'position', 'state']
        self.data[name] = None
    state = qbitproperty('state', 'Quantum state of qubit as 2-column')
    label = qbitproperty('label', 'Integer label asigned to qubit')
    position = qbitproperty('position', 'XYZ-coordinates')
    tag = qbitproperty('tag', 'Integer tag')
    #momentum = qbitproperty('momentum', 'XYZ-momentum')
    #mass = qbitproperty('mass', 'Atomic mass')
    #magmom = qbitproperty('magmom', 'Initial magnetic moment')
    #charge = qbitproperty('charge', 'Initial atomic charge')
    x = xyzproperty(0)
    y = xyzproperty(1)
    z = xyzproperty(2)

    #a = abcproperty(0)
    #b = abcproperty(1)
    #c = abcproperty(2)
