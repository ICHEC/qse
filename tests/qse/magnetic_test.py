import numpy as np
import pytest

import qse


@pytest.mark.parametrize("hsize, N", [(8, 4), (3, 3)])
def test_get_basis_shape(hsize, N):
    """Check that get_basis outputs a hsize * N shape"""
    ibasis = qse.magnetic.get_basis(hsize, N)
    assert ibasis.shape == (hsize, N)


@pytest.mark.parametrize(
    "statevector, ibasis, results",
    [
        (
            np.array([1 / np.sqrt(2), 1 / np.sqrt(2)], dtype=complex),
            qse.magnetic.get_basis(2**1, 1),
            np.array([1, 0, 0], dtype=complex).T.real,
        ),
        (
            np.array([1 / 2, 1 / 2, 1 / 2, 1 / 2], dtype=complex),
            qse.magnetic.get_basis(2**2, 2),
            np.array([[1, 0, 0], [1, 0, 0]], dtype=complex),
        ),
        (
            np.array([1 / np.sqrt(2), 0, 0, 1 / np.sqrt(2)], dtype=complex),
            qse.magnetic.get_basis(2**2, 2),
            np.array(
                [
                    [0, 0, 0],
                    [0, 0, 0],
                ],
                dtype=complex,
            ),
        ),
    ],
)
def test_get_spins_on_simple_statevectors(statevector, ibasis, results):
    """Test output of get_spins on simple quantum states"""
    N = int(np.log2(statevector.shape[0]))
    spins = qse.magnetic.get_spins(statevector, ibasis, N)
    assert np.allclose(results, spins)


@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_spin_values_are_less_than_one(k):
    "Test that the absolute values of the computed spins is less than 1"
    real_part = np.random.rand(2**k)
    imag_part = np.random.rand(2**k)
    statevector = real_part + 1j * imag_part
    statevector /= np.linalg.norm(statevector)
    ibasis = qse.magnetic.get_basis(2**k, k)
    spins = qse.magnetic.get_spins(statevector, ibasis, k)
    assert np.all(np.abs(spins) <= 1)
