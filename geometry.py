"""
geometry.py -- Geometric descriptor functor G : Struct_onc -> Geom

Implements coarse geometric descriptors used in the GTCF-Onc paper:
  - Kabsch superposition + backbone RMSD (used for KRAS switch I/II, Theorem 1)
  - radius of gyration / compactness descriptor (used for G-quadruplex, Theorem 2)
  - simple inter-residue "bend angle" descriptor (tubulin-style, reusable for
    any two-lobe conformational comparison)

These are intentionally coarse, closed-form descriptors -- matching the
"coarse, global compactness" scope condition stated in the paper's Theorem 2
and Limitations section. They are NOT claimed to be the only valid choice of
G; Section 8 of the paper states explicitly that a locally-resolved G could
behave differently.
"""
from __future__ import annotations
import numpy as np


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Backbone RMSD between two Nx3 coordinate sets after optimal
    (Kabsch) superposition. P and Q must already be in 1:1 residue
    correspondence (same length, same order)."""
    if P.shape != Q.shape:
        raise ValueError(f"shape mismatch: {P.shape} vs {Q.shape}")
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    Pc_rot = (R @ Pc.T).T
    diff = Pc_rot - Qc
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))


def radius_of_gyration(coords: np.ndarray) -> float:
    """Coarse compactness descriptor: Rg of a point cloud (e.g. all
    heavy atoms, or all P/glycosidic atoms of a quadruplex)."""
    c = coords - coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(c ** 2, axis=1))))


def pairwise_max_diameter(coords: np.ndarray) -> float:
    """Coarse size descriptor: largest pairwise distance in the point
    cloud (overall quadruplex/domain diameter)."""
    diffs = coords[:, None, :] - coords[None, :, :]
    d2 = np.sum(diffs ** 2, axis=-1)
    return float(np.sqrt(d2.max()))


def inter_lobe_bend_angle(lobe1: np.ndarray, lobe2: np.ndarray, hinge: np.ndarray) -> float:
    """Tubulin-style two-lobe bend angle (degrees), generalizable to any
    two-domain hinge geometry: angle between the vector from hinge to
    lobe1 centroid and hinge to lobe2 centroid."""
    v1 = lobe1.mean(axis=0) - hinge
    v2 = lobe2.mean(axis=0) - hinge
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_theta)))


def geometric_blindness_test(P: np.ndarray, Q: np.ndarray, threshold_angstrom: float = 1.0):
    """Definition 1/2 helper: returns (are_G_images_different, rmsd).
    G(S1) !~= G(S2) iff RMSD after optimal superposition exceeds threshold
    (default 1.0 A, well above the ~0.1-0.15 A coordinate error reported
    for high-resolution KRAS crystal structures such as 4OBE)."""
    rmsd = kabsch_rmsd(P, Q)
    return rmsd > threshold_angstrom, rmsd
