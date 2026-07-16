"""
topology.py -- Topological descriptor functor T : Struct_onc -> Top

Implements the persistence-diagram side of GTCF-Onc using the Vietoris-Rips
filtration (via the `ripser` package) on a Calpha or heavy-atom point cloud,
plus a bottleneck-distance comparison (via `persim`) used to operationalize
Definition 1 / Definition 2 (T(S1) ~= T(S2) vs T(S1) !~= T(S2)).

This is deliberately the same computational recipe as Section 6/7 of both
the base GTCF paper and the GTCF-Onc paper: Vietoris-Rips filtration,
Betti numbers (b0, b1, b2), bottleneck distance between diagrams.
"""
from __future__ import annotations
import numpy as np

try:
    from ripser import ripser
    from persim import bottleneck
    _HAVE_TDA = True
except ImportError:  # pragma: no cover
    _HAVE_TDA = False


def persistence_diagrams(coords: np.ndarray, maxdim: int = 2):
    """Compute persistence diagrams (H0, H1, H2) for a point cloud via a
    Vietoris-Rips filtration. Returns the list of diagrams as returned by
    ripser (one per homological dimension 0..maxdim)."""
    if not _HAVE_TDA:
        raise RuntimeError("ripser/persim not installed; pip install ripser persim")
    result = ripser(coords, maxdim=maxdim)
    return result["dgms"]


def betti_numbers(dgms, epsilon: float | None = None):
    """Count features alive at filtration value `epsilon` for each homology
    dimension. If epsilon is None, uses the median death value of H0 as a
    representative scale (a reasonable default for a single connected
    biomolecule point cloud)."""
    if epsilon is None:
        finite_deaths = dgms[0][:, 1]
        finite_deaths = finite_deaths[np.isfinite(finite_deaths)]
        epsilon = float(np.median(finite_deaths)) if len(finite_deaths) else 0.0
    betti = []
    for dgm in dgms:
        alive = np.sum((dgm[:, 0] <= epsilon) & (dgm[:, 1] > epsilon))
        betti.append(int(alive))
    return betti, epsilon


def bottleneck_distance(dgm_a, dgm_b):
    """Bottleneck distance between two persistence diagrams of the same
    homological dimension (Section 6/7 comparison metric)."""
    if not _HAVE_TDA:
        raise RuntimeError("ripser/persim not installed; pip install ripser persim")
    return float(bottleneck(dgm_a, dgm_b))


def topological_blindness_test(coords_a: np.ndarray, coords_b: np.ndarray,
                                maxdim: int = 2, threshold: float = 0.5):
    """Definition 1 helper: returns (are_T_images_equal, per-dimension
    bottleneck distances). T(S1) ~= T(S2) iff bottleneck distance in every
    tested dimension is below `threshold` (in the same length units as the
    input coordinates, e.g. Angstroms)."""
    dgms_a = persistence_diagrams(coords_a, maxdim=maxdim)
    dgms_b = persistence_diagrams(coords_b, maxdim=maxdim)
    distances = [bottleneck_distance(dgms_a[d], dgms_b[d]) for d in range(maxdim + 1)]
    are_equal = all(d < threshold for d in distances)
    return are_equal, distances
