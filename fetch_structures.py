"""
fetch_structures.py -- download the four witness structures from RCSB PDB.

NOTE ON PROVENANCE: this script could NOT be run to completion inside the
sandboxed environment that built this repository -- that environment's
network egress is restricted to a small allow-list of code-hosting domains
and does not include files.rcsb.org (confirmed: RCSB returns HTTP 403 to
requests from that sandbox). Run this script yourself, in any normal
internet-connected environment (a laptop, CI runner, Claude Code with
network access, etc.), to complete the retrieval and unblock
verify_theorem1.py / verify_theorem2.py end to end.

Usage:
    python src/fetch_structures.py
"""
from __future__ import annotations
import pathlib
import urllib.request

STRUCTURES = {
    "4OBE": "KRAS GDP-bound wild type (Theorem 1, S1)",
    "5VQ6": "KRAS GTP-gamma-S-bound wild type (Theorem 1, S2)",
    "143D": "human telomeric G-quadruplex, antiparallel basket, Na+ (Theorem 2, S3)",
    "1KF1": "human telomeric G-quadruplex, parallel propeller, K+ (Theorem 2, S4)",
}

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"


def fetch(pdb_id: str, out_dir: pathlib.Path) -> pathlib.Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pdb_id}.pdb"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Fetching {pdb_id} ({STRUCTURES[pdb_id]}) ...")
    urllib.request.urlretrieve(url, out_path)
    print(f"  -> saved to {out_path} ({out_path.stat().st_size} bytes)")
    return out_path


def main():
    for pdb_id in STRUCTURES:
        try:
            fetch(pdb_id, DATA_DIR / pdb_id)
        except Exception as exc:
            print(f"  FAILED for {pdb_id}: {exc}")
    print("\nDone. Re-run verify_theorem1.py and verify_theorem2.py.")


if __name__ == "__main__":
    main()
