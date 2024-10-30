"""
This module Provides QSE calculator interface to MyQLM. It is derived from
the ASE's calculator, though at the end it may look very different from
ASE calculator.
https://myqlm.github.io/
"""

import numpy as np
from qse.calc.calculator import Calculator, all_changes

qat_available = False
qlmaas_available = False

try:
    import qat.qpus
    if hasattr(qat.qpus, 'AnalogQPU'):
        AQPU_local = qat.qpus.AnalogQPU
        qat_available = True
    else:
        AQPU_local = None
except ModuleNotFoundError:
    qat_available = False

try:
    import qlmaas.qpus
    if hasattr(qlmaas.qpus, 'AnalogQPU'):
        AQPU_remote = qlmaas.qpus.AnalogQPU
        qlmaas_available = True
    else:
        AQPU_remote = None
except ModuleNotFoundError:
    qlmaas_available = False

if qat_available:
    AQPU = AQPU_local
else:
    if qlmaas_available:
        AQPU = AQPU_remote
    else:
        AQPU = None

if qat_available:
    print('qat is available')
else:
    print('qat is not available')
if qlmaas_available:
    print('qlmaas is available')
else:
    print('qlmaas is not available')
print(AQPU)

# analogQPU imported based on what's available

from qat.core.variables import Variable, heaviside

default_params = {'rydberg': {
        'amplitude': None,
        'detuning': None,
        'C6': 5420158.53, # unit (rad/µs)(µm)**6
        'min_omega': None,
        'max_omega': 12.57, # rad/µs from pulser
        'min_delta': None,
        'max_delta': 125.7, # rad/µs from pulser
        'min_atom_distance': 4, # µm from pulser
        'max_duration': 4000 # ns from pulser (max sequence duration)
        },
      'ssh': {},
      'fermion': {}
      }

# for rydberg system we need additional parameters: mix,max of amplitude and detuning to set based on device.
class Myqlm(Calculator):
    """QSE-Calculator for MyQLM
    """
    implemented_properties = ['energy', 'state', 'fidality']
    default_parameters = dict(label='q', qbits=None)

    def __init__(self, qbits=None,
                amplitude=None, detuning=None,
                qpu=None, analog=True,
                system='rydberg',
                label='myqlm-run'):
        super().__init__(label=label, qbits=qbits)
        self.qpu = AQPU if qpu is None else qpu
        self.label = label
        self.results = None
        self.system = system
        self.params = dict(default_params[self.system])
        # if any of the default parameters are passed
        # as arguments, use them to update self.params
        self.C6 = self.params['C6']
        #
    
    def occ_op(nqbits, qi):
        ti = qat.core.Term(1.0, 'Z', [qi])
        return (1 + qat.core.Observable(nqbits, pauli_terms=[ti])) / 2

    def _get_hamiltonian(self):
        if self.system == 'rydberg':
            ham = _generate_rydberg_hamiltonian()
        else:
            ham = None
        return ham
    #

    def _generate_rydberg_hamiltonian(self):
        rij = self.qbits.get_all_distances()
        nqbits = self.nqbits
        #> add checks for ensuring whether amplitudes/detunings etc
        #> are within allowed limits compatible with pulser virtual device
        amplitude = self.amplitude
        detuning = self.detuning

        H1_terms = [qat.core.Term(0.5, "X", [i]) for i in range(nqbits)]
        H2_terms = [qat.core.Term(0.5, "Z", [i]) for i in range(nqbits)]
        H_amplitude = qat.core.Observable(nqbits, pauli_terms=H1_terms)
        H_detuning = qat.core.Observable(nqbits, pauli_terms=H2_terms)
        H_interact = 0
        for i in range(nqbits):
            for j in range(i+1, nqbits):
                H_interact += (self.C6 / rij[i, j]**6) * occ_op(nqbits, i) * occ_op(nqbits, j)
        Hamiltonian = qat.core.Observable([
            (amplitude, H_amplitude),
            (detuning, H_detuning),
            (1, H_interact)
        ])
        return Hamiltonian
    #


    def calculate(self, qbits=None, properties=..., system_changes=...):
        #return super().calculate(qbits, properties, system_changes)
        pass


    


