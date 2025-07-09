import numpy as np
import pytest

import qse


@pytest.mark.parametrize(
    "hsize, N",
    [
        (8, 4),
        (3, 3)
    ]
)
def test_get_basis_returns_correct_qubit_space_shape(
    hsize, N
):
    ibasis = qse.magnetic.get_basis(hsize, N)
    assert ibasis.shape == (hsize, N)