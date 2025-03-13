import pytest
import qse
from qse.utils import squarelattice


def test_lattice():
    square_reg = squarelattice(1, 6, 6)
    assert isinstance(square_reg, qse.Qbits)


def test_lattice_fail():
    with pytest.raises(Exception):
        squarelattice('a', 2, 1)
