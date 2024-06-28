"""
This module defines abstract helper classes with the objective of reducing
boilerplace method definitions (i.e. duplication) in calculators.
"""

from abc import ABC, abstractmethod
from typing import Mapping, Any


class GetPropertiesMixin(ABC):
    """Mixin class which provides get_forces(), get_stress() and so on.

    Inheriting class must implement get_property()."""

    @abstractmethod
    def get_property(self, name, qbits=None, allow_calculation=True):
        """Get the named property."""

    #def get_energies(self, qbits=None):
    #    return self.get_property('energies', qbits)

    #def get_moment(self, qbits=None):
    #    return self.get_property('magmom', qbits)


""" class GetOutputsMixin(ABC):
    #Mixin class for providing get_fermi_level() and others.
    #
    #Effectively this class expresses data in calc.results as
    #methods such as get_fermi_level().
    #
    #Inheriting class must implement _outputmixin_get_results(),
    #typically returning self.results, which must be a mapping
    #using the naming defined in ase.outputs.Properties.
    
    @abstractmethod
    def _outputmixin_get_results(self) -> Mapping[str, Any]:
        #Return Mapping of names to result value.
        #This may be called many times and should hence not be
        #expensive (except possibly the first time).

    def _get(self, name):
        # Cyclic import, should restructure.
        from ase.calculators.calculator import PropertyNotPresent
        dct = self._outputmixin_get_results()
        try:
            return dct[name]
        except KeyError:
            raise PropertyNotPresent(name)

    def get_fermi_level(self):
        return self._get('fermi_level')

    def get_eigenvalues(self, kpt=0, spin=0):
        eigs = self._get('eigenvalues')
        return eigs[kpt, spin]

    def _eigshape(self):
        # We don't need this if we already have a Properties object.
        return self._get('eigenvalues').shape

    def get_occupation_numbers(self, kpt=0, spin=0):
        occs = self._get('occupations')
        return occs[kpt, spin]
 """
