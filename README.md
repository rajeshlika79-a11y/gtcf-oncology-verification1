# GTCF-Onc Computational Verification

Computational verification pipeline for the two oncology witness theorems in
*"Geometry–Topology Coequal Framework in Oncology (GTCF-Onc): Locating
Structural-Blindness and Topological-Sensitivity Gaps in Cancer Structural
Biology"* (Rajesh Kumar, independent researcher).

This repository implements Section 7 ("Required Computational Verification")
of that paper: it computes the geometric functor **G** and topological
functor **T** for the two witness pairs and tests whether they exhibit
Structural Blindness (Theorem 1) and Topological Sensitivity (Theorem 2).

| Theorem | Witness pair | Claim |
|---|---|---|
| 1 — Structural Blindness | KRAS GDP-bound (4OBE) vs KRAS GTP-γS-bound (5VQ6) | G differs (switch I/II rearrangement); T predicted equal |
| 2 — Topological Sensitivity | Telomeric G-quadruplex antiparallel basket (143D) vs parallel propeller (1KF1) | G equal (coarse compactness); T differs (loop connectivity) |

## Status of this repository, honestly stated

This repo was built inside a sandboxed agent environment whose network
egress is restricted to a small allow-list of code-hosting domains
(`github.com`, `pypi.org`, `npmjs.com`, etc.) and does **not** include
`files.rcsb.org`. That was verified directly: a request to
`files.rcsb.org` from that sandbox returns `HTTP 403 Forbidden`.

As a result:

- **The full pipeline code is complete, tested, and correct** — see
  `tests/test_pipeline.py` (7/7 passing) and the two smoke-test runs
  described in [`RESULTS.md`](RESULTS.md).
- **One real, partial fragment of PDB entry 4OBE** (header, full sequence,
  full secondary-structure records for both chains, and chain-A backbone
  atoms for residues 1–32) was retrieved through a web-fetch tool and is
  committed at `data/4OBE/4OBE_partial.pdb`. It is real crystallographic
  data, not synthetic — but it is one nucleotide state only, and it stops
  mid-structure (before switch II), so it cannot by itself establish or
  refute Theorem 1.
- **No coordinate data for 5VQ6, 143D, or 1KF1 was retrievable** in that
  sandboxed session.
- **Theorem 1 and Theorem 2 therefore remain conditional**, exactly as
  stated in the paper — this repository does not claim otherwise.

## How to complete the verification (takes ~1 minute with normal internet)

```bash
pip install -r requirements.txt
python src/fetch_structures.py      # downloads 4OBE, 5VQ6, 143D, 1KF1 from RCSB
python src/verify_theorem1.py       # full real-data run once files are present
python src/verify_theorem2.py       # full real-data run once files are present
```

Each `verify_theorem*.py` script auto-detects whether the full structure
files are present under `data/<PDB_ID>/<PDB_ID>.pdb`. If they are, it runs
the real Section 7 computation (Kabsch RMSD + Vietoris–Rips persistence
diagrams + bottleneck distance) and prints ESTABLISHED / NOT SUPPORTED for
that theorem. If they are not, it falls back to reporting what can
honestly be said from the partial real data plus a clearly-labeled
synthetic smoke test that only demonstrates the code is working correctly
— never presented as evidence for the paper's claims.

## Repository layout

```
data/                   PDB coordinate files (only 4OBE_partial.pdb is committed;
                        run fetch_structures.py to populate the rest)
src/
  pdb_utils.py          minimal PDB ATOM-record parser
  geometry.py           G functor: Kabsch RMSD, radius of gyration, diameter, bend angle
  topology.py           T functor: Vietoris-Rips persistence diagrams, bottleneck distance
  fetch_structures.py   downloads the four witness structures from RCSB PDB
  verify_theorem1.py    Structural Blindness verification (KRAS)
  verify_theorem2.py    Topological Sensitivity verification (G-quadruplex)
tests/
  test_pipeline.py      unit tests against synthetic ground truth (no network needed)
RESULTS.md              exact output of every script run so far, with dates
requirements.txt
LICENSE
```

## Citation

If you use this code, please cite the paper it accompanies:

> Kumar, R. *Geometry–Topology Coequal Framework in Oncology (GTCF-Onc):
> Locating Structural-Blindness and Topological-Sensitivity Gaps in Cancer
> Structural Biology.* Preprint, 2026.

## License

MIT (see `LICENSE`). PDB data retrieved via `fetch_structures.py` is
subject to RCSB PDB's own data usage policies.
