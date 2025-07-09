import numpy as np
import pytest

import qse

"""
@pytest.mark.parametrize(
    "hsize, N", [(8, 4), (3, 3)])
def test_get_basis_returns_correct_qubit_space_shape(hsize, N):
    ibasis = qse.magnetic.get_basis(hsize, N)
    assert ibasis.shape == (hsize, N)

"""

def add(x,y):
    print("x = ", x)
    print("y = ", y)
    print('-----------')

    return x + y 

@pytest.mark.parametrize(
    "x, y", [(2,2), (4,5)]
)
def test_add(x, y):
    assert add(x, y) == x ** 2 