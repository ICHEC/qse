import numpy as np
import pytest

import qse


def _qbits_checker(qbits, positions, total_qubits):
    assert isinstance(qbits, qse.Qbits)
    assert np.allclose(qbits.positions, positions)
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


def test_add_qbit():
    qbits = qse.Qbits(positions=np.arange(9).reshape(-1, 3), labels=[1, 2, 3])
    qbit = qse.Qbit(position=np.array([0.0, 0.0, 0.0]), label="q")
    qbits += qbit
    assert qbits.nqbits == 4
    assert all(qbits.labels == [1, 2, 3, "q"])


@pytest.mark.parametrize("nqbits", [2, 3, 4])
def test_get_distance(nqbits):
    """Test get_distance on a set of qbits."""
    positions = np.random.rand(nqbits, 3)
    qbits = qse.Qbits(positions=positions)

    assert qbits.nqbits == positions.shape[0]

    d = np.sqrt(((positions[0] - positions[1]) ** 2).sum())
    assert np.isclose(qbits.get_distance(0, 1), d)


@pytest.mark.parametrize("nqbits", [2, 3, 4])
def test_get_distances(nqbits):
    """Test get_distances on a set of qbits."""
    positions = np.random.rand(nqbits, 3)
    qbits = qse.Qbits(positions=positions)

    assert qbits.nqbits == positions.shape[0]
    inds = list(range(nqbits))
    d = np.sqrt(((positions[0] - positions) ** 2).sum(1))
    assert np.isclose(d[0], 0.0)
    assert np.allclose(qbits.get_distances(0, inds), d)


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


@pytest.mark.parametrize("nqbits", [2, 3, 4])
def test_get_all_distances_properties(nqbits):
    """Test the properties of get_all_distances."""
    qbits = qse.Qbits(positions=np.random.rand(nqbits, 3))
    distances = qbits.get_all_distances()
    assert distances.shape == (nqbits, nqbits)
    assert all(0.0 == i for i in np.diag(distances))
    assert np.allclose(distances, distances.T)


@pytest.mark.parametrize("nqbits", [3, 4])
@pytest.mark.parametrize("distance", [0.2, 2.1])
@pytest.mark.parametrize("inds", [[0, 1], [1, 2], [2, 0]])
def test_set_distance(nqbits, distance, inds):
    """Test set_distance."""
    qbits = qse.Qbits(positions=np.random.rand(nqbits, 3))
    qbits.set_distance(*inds, distance)
    assert np.isclose(qbits.get_distance(*inds), distance)


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

    assert qbits.positions.shape == positions.shape
    assert not np.allclose(qbits.positions, positions)


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
        assert np.allclose(qbits.positions, positions + disp)


def test_get_item():
    """Test __getitem__"""
    positions = np.random.rand(4, 3)
    qbits = qse.Qbits(positions=positions)

    # test int
    for indices in [0, 2]:
        assert isinstance(qbits[indices], qse.Qbit)
        assert np.allclose(qbits[indices].position, positions[indices])

    # test list
    for indices in [[0, 2], [1, 3, 2]]:
        assert isinstance(qbits[indices], qse.Qbits)
        assert np.allclose(qbits[indices].positions, positions[indices])

    # test slice
    assert isinstance(qbits[1:3], qse.Qbits)
    assert np.allclose(qbits[1:3].positions, positions[1:3])


@pytest.mark.parametrize("indices", [1, [0, 1, 3]])
def test_del_item(indices):
    """Test __delitem__"""
    nqbits = 4
    positions = np.random.rand(nqbits, 3)
    qbits = qse.Qbits(positions=positions)

    del qbits[indices]

    if isinstance(indices, int):
        indices = [indices]

    _qbits_checker(
        qbits,
        positions[[i for i in range(nqbits) if i not in indices]],
        nqbits - len(indices),
    )


def test_rotate():
    """Simple checks for rotate."""
    unit_square = np.array(
        [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
    )
    qbits = qse.Qbits(positions=unit_square)
    qbits.rotate(90, "z")

    assert not np.allclose(qbits.positions, unit_square)

    unit_square_rt = np.array(
        [[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    )
    assert np.allclose(qbits.positions, unit_square_rt)

    # Check that a centered square is invariant (with relabelling).
    square_centered = np.array(
        [[1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [-1.0, 1.0, 0.0], [-1.0, -1.0, 0.0]]
    )
    qbits = qse.Qbits(positions=square_centered)
    qbits.rotate(90, "z")
    assert np.allclose(qbits.positions, square_centered[[2, 0, 3, 1]])


@pytest.mark.parametrize("angle", [10, 20, 30])
@pytest.mark.parametrize("axis", [None, "z", (0, 0, 1)])
def test_rotate_square_z(angle, axis):
    """Test rotating a square system about the z axis."""
    positions = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0]]
    )
    qbits = qse.Qbits(positions=positions)

    if axis is None:
        qbits.rotate(angle)
    else:
        qbits.rotate(angle, axis)

    angle_rads = np.pi * angle / 180
    new_positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [np.cos(angle_rads), np.sin(angle_rads), 0.0],
            [0.0, 0.0, 1.0],
            [np.cos(angle_rads), np.sin(angle_rads), 1.0],
        ]
    )
    assert np.allclose(qbits.positions, new_positions)


@pytest.mark.parametrize("a", ["x", "-y", (0.0, 2.0, 3)])
@pytest.mark.parametrize("v", ["z", (1.0, 1.0, 0.0)])
@pytest.mark.parametrize("center", [(0.0, 0.0, 0.0), (-1, 3, 0.2)])
def test_rotate_distances(a, v, center):
    """Check a random rotation preserves distances."""
    positions = np.random.rand(4, 3)
    qbits = qse.Qbits(positions=positions)
    distances = qbits.get_all_distances()

    qbits.rotate(a, v, center)

    assert not np.allclose(qbits.positions, positions)
    assert np.allclose(qbits.get_all_distances(), distances)


def test_euler_rotate_and_rotate():
    """Test rotate and euler_rotate agree."""
    # Note that they rotate in different directions (clockwise & anti.)
    # May want to fix this.
    positions = np.random.rand(4, 3)

    qbits_1 = qse.Qbits(positions=positions)
    qbits_1.rotate(34, "z")

    qbits_2 = qse.Qbits(positions=positions)
    qbits_2.euler_rotate(-34)

    assert np.allclose(qbits_1.positions, qbits_2.positions)


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

    assert not np.allclose(qbits.positions, positions)
    assert np.allclose(qbits.get_all_distances(), distances)


@pytest.mark.parametrize("angle", [0.0, -44, -15.99, 5, 10, 36.9, 90, 180.0, 299.0])
def test_get_angle(angle):
    """Test get_angle on a simple 3-qbit system."""
    angle_rads = np.pi * angle / 180
    positions = np.array(
        [
            [np.cos(angle_rads), np.sin(angle_rads), 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
    )
    qbits = qse.Qbits(positions=positions)
    qbits.translate(10 * np.random.rand())  # random displacement
    qbits.euler_rotate(*(360 * np.random.rand(3)))  # random global rotation

    if angle > 180.0:
        angle = 360.0 - angle
    elif angle < 0.0:
        angle = -angle
    if np.isclose(angle, 0.0):
        # Tests often fail close to zero.
        assert np.isclose(qbits.get_angle(0, 1, 2), angle, atol=1e-5)
    else:
        assert np.isclose(qbits.get_angle(0, 1, 2), angle)


def test_get_angles():
    """Test get_angles."""
    indices = np.array(
        [
            [0, 1, 2],
            [0, 1, 3],
            [0, 1, 4],
            [0, 1, 2],
        ]
    )

    qbits = qse.Qbits(positions=np.random.rand(6, 3))
    angles = qbits.get_angles(indices)

    assert angles.shape == (indices.shape[0],)
    assert np.isclose(angles[0], angles[-1])


@pytest.mark.parametrize("angle", [11.2, 36.0])
@pytest.mark.parametrize("indices", [[0, 1, 2], [2, 4, 1], [0, 4, 2]])
def test_set_angle(angle, indices):
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
    assert not np.isclose(angle, qbits.get_angle(*indices))
    qbits.set_angle(*indices, angle)
    assert np.isclose(angle, qbits.get_angle(*indices))


@pytest.mark.parametrize("nqbits", [1, 2, 5])
@pytest.mark.parametrize(
    "repeats", [(1, 1, 1), (1, 2, 1), (1, 1, 3), (2, 3, 1), (2, 2, 2)]
)
def test_repeat(nqbits, repeats):
    """Test repeat (also testing imul / mul since they are all related.)."""
    cell = np.random.uniform(-2, 2, (3, 3))
    positions = np.random.uniform(-2, 2, (nqbits, 3))
    qbits = qse.Qbits(positions=positions, cell=cell)
    qbits = qbits.repeat(repeats)

    assert qbits.nqbits == nqbits * np.prod(repeats)

    ps = []
    for i in range(repeats[0]):
        for j in range(repeats[1]):
            for k in range(repeats[2]):
                v = i * cell[0] + j * cell[1] + k * cell[2]
                for p in positions:
                    print(p)
                    ps.append(v + p)

    ps = np.array(ps)
    assert np.allclose(ps, qbits.positions)


@pytest.mark.parametrize(
    "cell, expected",
    [
        (4.232 * np.eye(3), np.eye(3) / 4.232),  # square
        (
            np.array(  # rectangular
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 2.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0 / 2.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
        ),
        (
            np.array(  # hex
                [
                    [1.0, 0.0, 0.0],
                    [0.5, 0.5 * np.sqrt(3), 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            np.array(
                [
                    [1.0, -1.0 / np.sqrt(3), 0.0],
                    [0.0, 2.0 / np.sqrt(3), 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
        ),
        (
            0.5
            * np.array(  # fcc
                [
                    [0.0, 1.0, 1.0],
                    [1.0, 0.0, 1.0],
                    [1.0, 1.0, 0.0],
                ]
            ),
            np.array(
                [
                    [-1.0, 1.0, 1.0],
                    [1.0, -1.0, 1.0],
                    [1.0, 1.0, -1.0],
                ]
            ),
        ),
    ],
)
def test_reciprocal_cell(cell, expected):
    """Test the get_reciprocal_cell method."""
    qbits = qse.Qbits(cell=cell)
    r_cell = qbits.get_reciprocal_cell()
    assert np.allclose(r_cell, expected)

    # Check the two cells are orthonormal
    assert np.allclose(cell @ r_cell.T, np.eye(3))
