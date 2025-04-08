import numpy as np
import pytest

import qse


def _qbits_checker(qbits, positions, total_qubits):
    assert isinstance(qbits, qse.Qbits)
    assert np.allclose(qbits.get_positions(), positions)
    assert qbits.nqbits == total_qubits


@pytest.mark.parametrize("nqbits", [1, 2, 3, 10])
def test_basic_properties(nqbits):
    """Test basic Qbits methods."""
    positions = np.random.rand(nqbits, 3)
    qbits = qse.Qbits(positions=positions)
    _qbits_checker(qbits, positions, nqbits)
    assert len(qbits) == nqbits


@pytest.mark.parametrize("nqbits_1", [1, 2, 3])
@pytest.mark.parametrize("nqbits_2", [1, 2, 3])
def test_add(nqbits_1, nqbits_2):
    """Test adding Qbits together."""
    positions_1 = np.random.rand(nqbits_1, 3)
    qbits_1 = qse.Qbits(positions=positions_1)

    positions_2 = np.random.rand(nqbits_2, 3)
    qbits_2 = qse.Qbits(positions=positions_2)

    qbits_12 = qbits_1 + qbits_2

    _qbits_checker(
        qbits_12, np.concatenate((positions_1, positions_2)), nqbits_1 + nqbits_2
    )


def test_get_all_distances():
    """Test get_all_distances on a simple set of qbits."""
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
    )

    qbits = qse.Qbits(positions=positions)

    assert qbits.nqbits == positions.shape[0]

    distances = np.array(
        [
            [0.0, 1.0, 1.0],
            [1.0, 0.0, np.sqrt(2)],
            [1.0, np.sqrt(2), 0.0],
        ]
    )

    assert np.allclose(qbits.get_all_distances(), distances)


@pytest.mark.parametrize(
    "positions",
    [
        np.array([[-1, 0, 0], [1, 0, 0]]),
        np.array([[-1, 0, 0], [1, 0, 0], [1, 2, 0], [1, -2, 6], [-2, 0, -6]]),
    ],
)
@pytest.mark.parametrize("centroid", [0.0, 1.0, -1.0])
def test_get_centroid(positions, centroid):
    """Test get_centroid, note the positions parametrized all have zero centroid."""
    qbits = qse.Qbits(positions=positions + centroid)
    assert np.allclose(qbits.get_centroid(), [centroid] * 3)


@pytest.mark.parametrize(
    "positions",
    [
        np.array([[-1, 0, 0], [1, 0, 0]]),
        np.array([[-1, 0, 0], [1, 0, 0], [1, 2, 0], [1, -2, 6], [-2, 0, -6]]),
    ],
)
@pytest.mark.parametrize(
    "centroid", [[0.0] * 3, [1.0] * 3, [-1.0] * 3, np.random.rand(3)]
)
def test_set_centroid(positions, centroid):
    """Test set_centroid, note the positions parametrized all have zero centroid."""
    qbits = qse.Qbits(positions=positions)
    qbits.set_centroid(centroid)
    assert np.allclose(qbits.get_centroid(), centroid)


@pytest.mark.parametrize("nqbits", [1, 2, 3, 10])
def test_rattle(nqbits):
    """Test rattle."""
    positions = np.random.rand(nqbits, 3)
    qbits = qse.Qbits(positions=positions)
    qbits.rattle()

    assert qbits.get_positions().shape == positions.shape
    assert not np.allclose(qbits.get_positions(), positions)


@pytest.mark.parametrize("nqbits", [1, 2, 3, 10])
@pytest.mark.parametrize(
    "type_of_disp", ["scalar", "vector", "matrix", "vector_error", "matrix_error"]
)
def test_translate(nqbits, type_of_disp):
    """Test translate."""
    positions = np.random.rand(nqbits, 3)
    qbits = qse.Qbits(positions=positions)

    shape_dict = {
        "scalar": (),
        "vector": (3,),
        "matrix": (nqbits, 3),
        "vector_error": (4,),
        "matrix_error": (nqbits + 1, 3),
    }
    disp = np.random.rand(*shape_dict[type_of_disp])

    if "error" in type_of_disp:
        with pytest.raises(Exception):
            qbits.translate(disp)
    else:
        qbits.translate(disp)
        assert np.allclose(qbits.get_positions(), positions + disp)


def test_get_item():
    """Test __getitem__"""
    positions = np.random.rand(4, 3)
    qbits = qse.Qbits(positions=positions)

    # test int
    for indicies in [0, 2]:
        assert isinstance(qbits[indicies], qse.Qbit)
        assert np.allclose(qbits[indicies].position, positions[indicies])

    # test list
    for indicies in [[0, 2], [1, 3, 2]]:
        assert isinstance(qbits[indicies], qse.Qbits)
        assert np.allclose(qbits[indicies].get_positions(), positions[indicies])

    # test slice
    assert isinstance(qbits[1:3], qse.Qbits)
    assert np.allclose(qbits[1:3].get_positions(), positions[1:3])


@pytest.mark.parametrize("indicies", [1, [0, 1, 3]])
def test_get_item(indicies):
    """Test __delitem__"""
    nqbits = 4
    positions = np.random.rand(nqbits, 3)
    qbits = qse.Qbits(positions=positions)

    del qbits[indicies]

    if isinstance(indicies, int):
        indicies = [indicies]

    _qbits_checker(
        qbits,
        positions[[i for i in range(nqbits) if i not in indicies]],
        nqbits - len(indicies),
    )


def test_rotate():
    """Simple checks for rotate."""
    unit_square = np.array(
        [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
    )
    qbits = qse.Qbits(positions=unit_square)
    qbits.rotate(90, "z")

    assert not np.allclose(qbits.get_positions(), unit_square)

    unit_square_rt = np.array(
        [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    )
    assert np.allclose(qbits.get_positions(), unit_square_rt)

    # Check that a centered square is invariant (with relabelling).
    square_centered = np.array(
        [[1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [-1.0, 1.0, 0.0], [-1.0, -1.0, 0.0]]
    )
    qbits = qse.Qbits(positions=square_centered)
    qbits.rotate(90, "z")
    assert np.allclose(qbits.get_positions(), square_centered[[2, 0, 3, 1]])


@pytest.mark.parametrize("angle", [10, 20, 30])
def test_rotation_square_z(angle):
    positions = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0]]
    )
    qbits = qse.Qbits(positions=positions)
    qbits.rotate(angle, "z")

    new_positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [np.cos(angle), np.sin(angle), 0.0],
            [0.0, 0.0, 1.0],
            [np.cos(angle), np.sin(angle), 1.0],
        ]
    )
    assert np.allclose(qbits.get_positions(), new_positions)


@pytest.mark.parametrize("a", ["x", "-y", (0.0, 2.0, 3)])
@pytest.mark.parametrize("v", ["z", (1.0, 1.0, 0.0)])
@pytest.mark.parametrize("center", [(0.0, 0.0, 0.0), (-1, 3, 0.2)])
def test_rotate_distances(a, v, center):
    """Check a random rotation preserves distances."""
    positions = np.random.rand(4, 3)
    qbits = qse.Qbits(positions=positions)
    distances = qbits.get_all_distances()

    qbits.rotate(a, v, center)

    assert not np.allclose(qbits.get_positions(), positions)
    assert np.allclose(qbits.get_all_distances(), distances)


def test_euler_rotate_and_rotate():
    """Test rotate and euler_rotate agree."""
    ### Note that they rotate in different directions (clockwise & anti.)
    ### May want to fix this.
    positions = np.random.rand(4, 3)

    qbits_1 = qse.Qbits(positions=positions)
    qbits_1.rotate(34, "z")

    qbits_2 = qse.Qbits(positions=positions)
    qbits_2.euler_rotate(-34)

    assert np.allclose(qbits_1.get_positions(), qbits_2.get_positions())


@pytest.mark.parametrize("phi", [11.2, 45.0])
@pytest.mark.parametrize("theta", [-3.0, 40])
@pytest.mark.parametrize("psi", [18.2])
@pytest.mark.parametrize("center", [(0.0, 0.0, 0.0), (-1, 3, 0.2)])
def test_euler_rotate_distances(phi, theta, psi, center):
    """Check a random euler rotation preserves distances."""
    positions = np.random.rand(4, 3)
    qbits = qse.Qbits(positions=positions)
    distances = qbits.get_all_distances()

    qbits.euler_rotate(phi, theta, psi, center)

    assert not np.allclose(qbits.get_positions(), positions)
    assert np.allclose(qbits.get_all_distances(), distances)


@pytest.mark.parametrize("angle", [90, 36.0])
def test_get_angle(angle):
    """Test get_angle on a simple 3-qbit system."""
    angle_rads = np.pi * angle / 180
    positions = np.array(
        [
            [np.cos(angle_rads), np.sin(angle_rads), 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    qbits = qse.Qbits(positions=positions)
    assert np.isclose(qbits.get_angle(0, 1, 2), angle)

    # Angle should be invariant under global translation.
    qbits.translate(-3.3)
    assert np.isclose(qbits.get_angle(0, 1, 2), angle)

    # Angle should be invariant under global rotation.
    qbits.euler_rotation(10, 20, -44.5)
    assert np.isclose(qbits.get_angle(0, 1, 2), angle)


@pytest.mark.parametrize("angle", [11.2, 36.0])
@pytest.mark.parametrize("indicies", [[0, 1, 2], [2, 4, 1], [0, 4, 2]])
def test_set_angle(angle, indicies):
    """Test set_angle on a simple 5-qbit system."""
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ]
    )
    qbits = qse.Qbits(positions=positions)
    assert not np.isclose(angle, qbits.get_angle(*indicies))
    qbits.set_angle(*indicies, angle)
    assert np.isclose(angle, qbits.get_angle(*indicies))
