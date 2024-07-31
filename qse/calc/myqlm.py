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
    qat_available = True
except ModuleNotFoundError:
    qat_available = False
try:
    import qlmaas.qpus
    qlmaas_available = True
except ModuleNotFoundError:
    qlmaas_available = False

if qat_available:
    if hasattr(qat.qpus, 'AnalogQPU'):
        AQPU = qat.qpus.AnalogQPU
    else:
        AQPU = None
        qat_available = False
if qlmaas_available and not qat_available:
    if hasattr(qlmaas.qpus, 'AnalogQPU'):
        AQPU = qlmaas.qpus.AnalogQPU
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


default_params = {'rydberg': {
        'amplitude': None,
        'detuning': None,
        'C6': None,
        'min_omega': None,
        'max_omega': None,
        'min_delta': None,
        'max_delta': None,
        'min_atom_distance': None
        },
      'ssh': {},
      'heisenberg': {}
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
        self.parameters = None
        self.results = None
        self.system = system
        default = default_params[self.system]
    

    def calculate(self, qbits=None, properties=..., system_changes=...):
        #return super().calculate(qbits, properties, system_changes)
        pass

    


