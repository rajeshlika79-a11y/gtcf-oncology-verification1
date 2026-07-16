# Results Log

This file records the exact, unedited output of every verification run
performed on this pipeline, with the environment and date noted. Update it
by appending — do not delete prior entries — whenever you (re)run the
scripts, especially after `fetch_structures.py` succeeds and the full-data
branch of `verify_theorem1.py` / `verify_theorem2.py` executes for the
first time.

---

## Run 1 — 2026-07-16, sandboxed build environment (no RCSB network access)

Environment: Python 3.12, numpy/scipy/ripser/persim/biopython installed via
pip. Network egress restricted to an allow-list that excludes
`files.rcsb.org` (confirmed via direct request -> `HTTP 403 Forbidden`).

### Unit tests

```
$ python -m pytest tests/ -v
tests/test_pipeline.py::test_kabsch_rmsd_zero_for_pure_rigid_motion PASSED
tests/test_pipeline.py::test_kabsch_rmsd_positive_for_internal_change PASSED
tests/test_pipeline.py::test_radius_of_gyration_scales_with_spread PASSED
tests/test_pipeline.py::test_pairwise_max_diameter_matches_known_pair PASSED
tests/test_pipeline.py::test_geometric_blindness_threshold_logic PASSED
tests/test_pipeline.py::test_persistence_diagram_detects_loop PASSED
tests/test_pipeline.py::test_bottleneck_distance_zero_for_identical_diagrams PASSED
7 passed in 1.38s
```

**Interpretation:** the geometry and topology code (Kabsch RMSD, radius of
gyration, pairwise diameter, Vietoris–Rips persistence diagrams, bottleneck
distance) is verified correct against synthetic point clouds with known
ground truth, independent of any PDB data availability.

### Theorem 1 (KRAS) — `python src/verify_theorem1.py`

Full `4OBE.pdb` / `5VQ6.pdb` were not present (network-blocked in this
environment). Fallback branch executed:

```
=== PARTIAL REAL DATA (Theorem 1 input data not fully retrieved) ===
Real partial data available: 4OBE chain A, residues 1-32 (32 residues, 32 resolved Calpha atoms).
  Radius of gyration of this fragment (P-loop + start of switch I): 11.08 A

--- SYNTHETIC SMOKE TEST (code-correctness check only; NOT real evidence) ---
  [synthetic] G-images RMSD = 2.047 A -> DIFFER
  [synthetic] T-image bottleneck distances = [0.093, 0.128] -> MATCH (as predicted by Theorem 1)
```

**Interpretation:** the ONLY real number here is the radius of gyration of
the actual 4OBE chain-A residues 1–32 fragment (11.08 Å). The RMSD/
bottleneck-distance numbers are from a synthetic hinge-rotation constructed
to mimic the order-of-magnitude switch-I/II rearrangement reported in the
literature; they demonstrate the code reproduces the qualitative pattern
Theorem 1 predicts (G differs, T unchanged), but are **not** a computation
on real 4OBE/5VQ6 coordinates and must not be cited as such.
**Theorem 1 status: still CONDITIONAL**, unchanged from the paper.

### Theorem 2 (G-quadruplex) — `python src/verify_theorem2.py`

No real 143D/1KF1 coordinate data was retrievable in this environment.
Fallback branch executed on a fully synthetic idealized tetrad-stack
construction (see `_idealized_tetrad_stack`, `_basket_topology_loops`,
`_propeller_topology_loops` in `src/verify_theorem2.py`):

```
=== SYNTHETIC SMOKE TEST (code-correctness check only; NOT real 143D/1KF1 data) ===
[synthetic] G(basket): Rg = 8.59, max diameter = 19.50
[synthetic] G(propeller): Rg = 9.45, max diameter = 21.03
[synthetic] G-images MATCH (coarse compactness)
[synthetic] T-image bottleneck distances = [3.118, 2.343] -> DIFFER (as predicted by Theorem 2)
```

**Interpretation:** entirely synthetic — no real 143D/1KF1 atoms were used
anywhere in this computation. It shows the code correctly distinguishes
loop-connectivity topology from coarse compactness when the two are
constructed to match/differ by design. **Theorem 2 status: still
CONDITIONAL**, unchanged from the paper.

---

## Run 2 — (fill in after running `fetch_structures.py` with real internet access)

```
$ python src/fetch_structures.py
...

$ python src/verify_theorem1.py
...

$ python src/verify_theorem2.py
...
```

Once this run is filled in with real full-file output, update the theorem
status lines above and in the paper (Section 7 / Section 8) accordingly:
if `verify_theorem1.py` reports `THEOREM 1 STATUS: ESTABLISHED` with a
bottleneck distance below threshold, Theorem 1 moves from conditional to
established, and the same for Theorem 2.
