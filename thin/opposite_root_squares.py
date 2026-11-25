"""Check the opposite-root cycle systems at one line of the Ree–Tits model."""

from collections import Counter

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path

from check_opposite_seed import ROOT, octagon


def radial_walks(target_is_line):
    """Counts of nonbacktracking incidence walks ending at a fixed target."""
    shells = np.arange(9)
    degree = np.where((shells + int(target_is_line)) % 2 == 0, 5, 3)
    inward = np.ones(9, dtype=int)
    inward[0], inward[8] = 0, degree[8]
    radial = np.diag(inward[1:], -1) + np.diag((degree - inward)[:-1], 1)
    at_target = (shells == 0).astype(int)
    walks = [at_target, radial @ at_target]
    walks.append(radial @ walks[1] - degree * at_target)
    for length in range(3, 11):
        walks.append(radial @ walks[length - 1] - (degree - 1) * walks[length - 2])
    return walks


def attachment_counts():
    walks = [radial_walks(False), radial_walks(True)]
    counts = []
    for point_distance in (6, 8):
        # Trim forbidden edges alternately at the start and end of the walk.
        distances = ((point_distance, 7), (7, 8))
        total = sum((-1)**(i+j) * walks[j % 2][10-i-j][distances[i % 2][j % 2]]
                    for i in range(11) for j in range(11-i))
        counts.append(int(total))
    assert counts == [18, 13]
    return np.where(np.eye(3, dtype=bool), counts[0], counts[1])


def component_lengths(distance, points, lines, x, y, expected_blocks):
    assert distance[x, y] == 4
    left = np.flatnonzero((distance[x] == 2) & (distance[y] == 4))
    right = np.flatnonzero((distance[y] == 2) & (distance[x] == 4))
    incidence = (distance[np.ix_(left, right)] == 2).astype(np.int16)
    assert incidence.shape == (66, 66)
    assert (incidence.sum(axis=0) == 2).all() and (incidence.sum(axis=1) == 2).all()

    x_points = lines[x]
    projections = points[np.ix_(x_points, lines[y])]
    assert ((projections == 3).sum(axis=1) == 1).all()
    assert ((projections == 3).sum(axis=0) == 1).all()
    y_points = lines[y][projections.argmin(axis=1)]
    left_color = points[x_points[:, None, None], lines[left]].min(axis=2).argmin(axis=0)
    right_color = points[y_points[:, None, None], lines[right]].min(axis=2).argmin(axis=0)
    rows, cols = np.nonzero(incidence)
    blocks = np.zeros((3, 3), dtype=int)
    np.add.at(blocks, (left_color[rows], right_color[cols]), 1)
    assert np.array_equal(blocks, expected_blocks)

    # Two steps stay on the left. Each cycle has equally many vertices on each side.
    _, component = connected_components(csr_matrix(incidence @ incidence.T), directed=False)
    return tuple(sorted(map(int, 2 * np.bincount(component))))


def main():
    expected_blocks = attachment_counts()
    points, lines = octagon(np.load(ROOT / "data/generators.npy"))
    incidence = csr_matrix(
        (np.ones(lines.size, dtype=int), (np.repeat(np.arange(len(lines)), 3), lines.ravel())),
        shape=(2925, 1755),
    )
    distance = shortest_path(incidence @ incidence.T, unweighted=True).astype(np.int8)
    types = Counter(component_lengths(distance, points, lines, 0, y, expected_blocks)
                    for y in np.flatnonzero(distance[0] == 4))
    assert types == {(6, 12, 12, 14, 14, 14, 30, 30): 2048}
    print("Cycle lengths and partner counts at line 0", dict(types))


if __name__ == "__main__":
    main()
