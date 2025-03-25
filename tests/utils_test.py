import numpy as np
import pytest

import qse


@pytest.mark.parametrize("lattice_spacing", [0.5, 3.1])
@pytest.mark.parametrize("N1", [2, 3])
@pytest.mark.parametrize("N2", [2, 3])
def test_squarelattice(lattice_spacing, N1, N2):
    qbits = qse.utils.squarelattice(lattice_spacing, N1, N2)
    assert isinstance(qbits, qse.Qbits)
    assert qbits.nqbits == N1 * N2

    # Check the distance between the first qubit and it's nearest neighbour
    # is the lattice spacing.
    assert np.isclose(qbits.get_all_distances()[0][1:].min(), lattice_spacing)
