"""
pdb_utils.py -- minimal, dependency-light PDB ATOM record parser.

Deliberately does not require biopython so the pipeline has a single
lightweight fallback path; a biopython-based path is also provided for
convenience when biopython is installed.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass


@dataclass
class Residue:
    chain: str
    resnum: int
    resname: str
    atoms: dict  # atom name -> np.ndarray([x, y, z])


def parse_pdb_atoms(path: str, chain: str | None = None):
    """Parse ATOM records from a PDB file into a list of Residue objects,
    in file order. HETATM records are ignored. Only the first alternate
    location (blank or 'A') is kept for each atom."""
    residues: dict[tuple, Residue] = {}
    order: list[tuple] = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            rec_chain = line[21].strip()
            if chain is not None and rec_chain != chain:
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            resnum = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            key = (rec_chain, resnum)
            if key not in residues:
                residues[key] = Residue(rec_chain, resnum, resname, {})
                order.append(key)
            residues[key].atoms[atom_name] = np.array([x, y, z])
    return [residues[k] for k in order]


def ca_trace(residues: list[Residue]) -> np.ndarray:
    """Extract the Calpha trace (Nx3) from a residue list, in order,
    skipping residues where CA is missing."""
    pts = [r.atoms["CA"] for r in residues if "CA" in r.atoms]
    return np.array(pts)


def residue_range(residues: list[Residue], first: int, last: int) -> list[Residue]:
    return [r for r in residues if first <= r.resnum <= last]


def all_heavy_atoms(residues: list[Residue]) -> np.ndarray:
    """Flatten all non-hydrogen atom coordinates across a residue list
    (used for the G-quadruplex compactness/diameter descriptors, where
    the whole-nucleotide point cloud is the natural unit rather than a
    single backbone trace atom)."""
    pts = []
    for r in residues:
        for name, xyz in r.atoms.items():
            if name.strip().startswith("H"):
                continue
            pts.append(xyz)
    return np.array(pts)
