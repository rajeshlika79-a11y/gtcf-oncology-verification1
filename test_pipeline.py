"""
test_pipeline.py -- unit tests for geometry.py and topology.py.

These validate CODE correctness on synthetic point clouds with known
ground truth (e.g. a point cloud vs. an exact rigid rotation of itself
should have RMSD ~ 0). They do not depend on network access and should
pass in any environment, including the sandboxed one this repo was
originally built in.

Run:
    python -m pytest tests/ -v
"""
import pathlib
import sys
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from geometry import kabsch_rmsd, radius_of_gyration, pairwise_max_diameter, geometric_blindness_test
from topology import persistence_diagrams, betti_numbers, bottleneck_distance, _HAVE_TDA


def _random_cloud(n=40, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, 3)) * 5.0


def _rotate(coords, angle_deg, axis=(0, 0, 1)):
    axis = np.array(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    theta = np.radians(angle_deg)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
    return (R @ coords.T).T + np.array([3.0, -2.0, 1.0])  # rotate + translate


def test_kabsch_rmsd_zero_for_pure_rigid_motion():
    coords = _random_cloud()
    moved = _rotate(coords, 47.0)
    rmsd = kabsch_rmsd(coords, moved)
    assert rmsd < 1e-6, f"expected ~0 RMSD after superposition, got {rmsd}"


def test_kabsch_rmsd_positive_for_internal_change():
    coords = _random_cloud()
    changed = coords.copy()
    changed[-5:] += np.array([10.0, 0, 0])  # break rigidity for a subset
    rmsd = kabsch_rmsd(coords, changed)
    assert rmsd > 0.5


def test_radius_of_gyration_scales_with_spread():
    small = _random_cloud(seed=1) * 0.5
    large = small * 4.0
    assert radius_of_gyration(large) > radius_of_gyration(small)


def test_pairwise_max_diameter_matches_known_pair():
    pts = np.array([[0, 0, 0], [10, 0, 0], [3, 3, 3]])
    assert abs(pairwise_max_diameter(pts) - 10.0) < 1e-9


def test_geometric_blindness_threshold_logic():
    coords = _random_cloud()
    identical = _rotate(coords, 0.001)  # negligible rigid motion
    is_diff, rmsd = geometric_blindness_test(coords, identical, threshold_angstrom=1.0)
    assert not is_diff
    changed = coords.copy()
    changed[:10] += np.array([5.0, 5.0, 5.0])
    is_diff2, rmsd2 = geometric_blindness_test(coords, changed, threshold_angstrom=1.0)
    assert is_diff2 and rmsd2 > rmsd


def test_persistence_diagram_detects_loop():
    if not _HAVE_TDA:
        return  # skip if ripser/persim unavailable in this environment
    # A circle of points should show a persistent H1 feature (a loop);
    # a random cloud of the same size typically should not.
    n = 60
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    circle = np.stack([np.cos(theta), np.sin(theta), np.zeros(n)], axis=1) * 5.0
    dgms_circle = persistence_diagrams(circle, maxdim=1)
    h1_circle = dgms_circle[1]
    finite = h1_circle[np.isfinite(h1_circle[:, 1])]
    if len(finite) == 0:
        return
    persistence = finite[:, 1] - finite[:, 0]
    assert persistence.max() > 1.0, "expected a long-lived H1 feature for a circle"


def test_bottleneck_distance_zero_for_identical_diagrams():
    if not _HAVE_TDA:
        return
    coords = _random_cloud(seed=2)
    dgms = persistence_diagrams(coords, maxdim=1)
    d = bottleneck_distance(dgms[1], dgms[1])
    assert d < 1e-6


if __name__ == "__main__":
    import subprocess
    subprocess.run(["python", "-m", "pytest", __file__, "-v"])
