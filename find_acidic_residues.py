"""
find_acidic_residues.py
=======================
Identifies the two catalytic acidic residues (Glu/Asp) nearest to a
carbohydrate ligand in a batch of PDB files.

Strategy
--------
1. Parse each PDB file with Biopython.
2. Locate all residues that are carbohydrate-like (via a curated
   residue-name list OR by atoms that are not common small molecules).
3. Compute the geometric centroid of the carbohydrate heavy atoms.
4. Collect all Glu and Asp residues; for each, take the representative
   "functional" atom (OE1/OE2 for Glu, OD1/OD2 for Asp).
5. Rank by minimum distance from any carb atom; return the two closest.
6. Optionally validate against a reference structure (superposition-free,
   centroid-anchored) to filter out false positives.

Requirements
------------
    pip install biopython numpy pandas tqdm

Usage
-----
    # Basic: process a folder of PDB files
    python find_acidic_residues.py --pdb_dir /path/to/pdbs --out results.csv

    # With reference structure to tune the search radius
    python find_acidic_residues.py --pdb_dir /path/to/pdbs \
        --reference ref.pdb --ref_carb_chain A --out results.csv

    # Explicit carbohydrate residue name(s) (comma-separated 3-letter codes)
    python find_acidic_residues.py --pdb_dir /path/to/pdbs \
        --carb_names MAN,FUC,GLA,BGC,GLC --out results.csv
"""

import argparse
import csv
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

# Suppress Biopython warnings about non-standard records
warnings.filterwarnings("ignore")

try:
    from Bio import PDB
    from Bio.PDB import PDBParser, NeighborSearch
except ImportError:
    sys.exit("Biopython is required:  pip install biopython")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ACIDIC_RESIDUES = {"GLU", "ASP"}

# Carboxylic / side-chain oxygen atoms used as the functional point
ACIDIC_OXY = {
    "GLU": ["OE1", "OE2"],
    "ASP": ["OD1", "OD2"],
}

# Common solvent / ion residue names to exclude from "carbohydrate" search
EXCLUDE_HETS = {
    "HOH", "WAT", "H2O",                          # water
    "SO4", "PO4", "CL", "NA", "MG", "ZN", "CA",   # ions/salts
    "EDO", "GOL", "PEG", "MPD", "FMT", "ACT",     # cryo-protectants
    "BME", "DTT", "TRS", "MES", "HEP", "EPE",     # buffers
    "ATP", "ADP", "AMP", "NAD", "FAD", "FMN",     # cofactors
    "HEM", "CLA", "BCL",                           # porphyrins
}

# Well-known monosaccharide / sugar residue names (extend as needed)
KNOWN_CARBS = {
    # Hexoses
    "GLC", "BGC", "BGLC", "GLA", "GAL", "MAN", "BMA",
    "FRU", "TAG", "SOR",
    # Pentoses
    "RIB", "ARA", "XYL", "LYX",
    # Deoxy sugars
    "FUC", "FCA", "RHA", "RHM",
    # GlcNAc / GalNAc / ManNAc
    "NAG", "NDG", "NGA", "A2G", "BNG", "MNA",
    # Neuraminic acids
    "SIA", "SLB", "NGC",
    # Disaccharides (sometimes appear as single residues)
    "LAT", "MAL", "CEL", "TRE",
    # Uronic acids
    "GCU", "IDR",
    # Generic / other
    "SGN", "SHB", "GCS", "l:b", "LSB", "LB"

}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_atoms(residue):
    """Return a list of atoms in a Biopython Residue object."""
    return list(residue.get_atoms())


def centroid(atoms):
    """Compute the 3-D centroid of a list of Biopython Atom objects."""
    coords = np.array([a.get_vector().get_array() for a in atoms])
    return coords.mean(axis=0)


def min_dist_to_set(atom, atom_set):
    """
    Minimum Euclidean distance from a single atom to a collection of atoms.
    Uses Biopython's subtraction operator (returns a Vector).
    """
    v = atom.get_vector()
    dists = [np.sqrt(np.sum((v - b.get_vector()), 2)) for b in atom_set]
    return min(dists) if dists else float("inf")


def is_carbohydrate(residue, carb_names):
    """
    Decide whether a residue is a carbohydrate.
    Comparison is always case-insensitive (names are uppercased before lookup).
    If carb_names is non-empty, use that whitelist only.
    Otherwise fall back to KNOWN_CARBS and exclude EXCLUDE_HETS.
    """
    resname = residue.get_resname().strip().upper()
    if carb_names:
        # Normalise to uppercase in case user typed lowercase (e.g. "l:b")
        return resname in {n.upper() for n in carb_names}
    if resname in EXCLUDE_HETS:
        return False
    return resname in KNOWN_CARBS



def collect_carb_atoms_raw(pdb_path, carb_names):
    """
    Parse ATOM/HETATM lines with fixed-width column slicing (PDB standard)
    instead of relying on Biopython, which misparsed residue names containing
    non-alphanumeric characters (e.g. "l:b").

    Returns a list of numpy arrays (xyz coordinates) — not Biopython Atom
    objects — so the caller must use raw numpy distance math instead of
    Biopython Vector arithmetic.

    PDB column layout (1-based, inclusive):
      cols  1- 6  record type  (ATOM / HETATM)
      cols  7-11  serial
      cols 13-16  atom name
      cols 17     alt loc
      cols 18-20  residue name   ← 3 chars
      col  21     residue name 4th char (non-standard; sometimes used)
      col  22     chain ID
      cols 23-26  residue seq number
      col  27     insertion code
      cols 31-38  x
      cols 39-46  y
      cols 47-54  z
    """
    carb_names_upper = {n.strip().upper() for n in carb_names}
    coords = []
    with open(pdb_path, "r", errors="replace") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            if len(line) < 54:
                continue
            # Residue name: cols 18-21 (0-based 17:21), strip and upper
            resname = line[17:21].strip().upper()
            atom_name = line[12:16].strip()
            # Skip hydrogen atoms
            if atom_name.startswith("H") or atom_name.startswith("1H") or atom_name.startswith("2H"):
                continue
            if resname in carb_names_upper:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append(np.array([x, y, z]))
                except ValueError:
                    continue
    return coords


def min_dist_raw(coord, carb_coords):
    """Minimum distance from a numpy xyz array to a list of numpy xyz arrays."""
    if not carb_coords:
        return float("inf")
    diffs = np.array(carb_coords) - coord
    return float(np.sqrt(np.sum(np.pow(diffs, 2))))

def collect_carb_atoms(structure, carb_names):
    """
    Return all heavy atoms that belong to carbohydrate residues.

    Accepts BOTH HETATM (het_flag starts with 'H_') AND plain ATOM records
    (het_flag == ' ') whose residue name matches the carbohydrate whitelist.
    Some force-field / design outputs (e.g. RFdiffusion) write ligands as
    ATOM records instead of HETATM, so filtering on het_flag alone misses them.
    """
    carb_atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                het_flag = res.get_id()[0]
                is_het  = het_flag.startswith("H_")
                is_atom = het_flag == " "
                if (is_het or is_atom) and is_carbohydrate(res, carb_names):
                    for atom in res.get_atoms():
                        if not atom.get_name().startswith("H"):   # skip H atoms
                            carb_atoms.append(atom)
    return carb_atoms


def collect_acidic_residues(structure):
    """
    Return all (residue_obj, chain_id, res_id) for Glu/Asp residues
    in ATOM records (standard amino acids).
    """
    acidic = []
    for model in structure:
        for chain in model:
            for res in chain:
                het_flag = res.get_id()[0]
                resname = res.get_resname().strip().upper()
                if het_flag == " " and resname in ACIDIC_RESIDUES:
                    acidic.append((res, chain.get_id(), res.get_id()))
    return acidic


def functional_atoms(residue):
    """
    Return the side-chain oxygen atoms of a Glu or Asp residue.
    Falls back to all atoms if named atoms are absent.
    """
    resname = residue.get_resname().strip().upper()
    target_names = ACIDIC_OXY.get(resname, [])
    atoms = [residue[n] for n in target_names if n in residue]
    if not atoms:
        atoms = list(residue.get_atoms())
    return atoms


def find_two_closest_acidic(structure, carb_atoms, cutoff=12.0):
    """
    Find the two acidic residues whose side-chain oxygens are closest
    to any carbohydrate atom, within `cutoff` Ångströms.

    Parameters
    ----------
    cutoff : float
        Initial search sphere radius in Å.  If fewer than 2 residues are
        found, the cutoff is expanded automatically (up to 30 Å).

    Returns
    -------
    list of dict  (length 0, 1, or 2)
    """
    acidic = collect_acidic_residues(structure)
    if not acidic:
        return []

    scored = []
    for res, chain_id, res_id in acidic:
        f_atoms = functional_atoms(res)
        d = min(min_dist_to_set(fa, carb_atoms) for fa in f_atoms)
        scored.append({
            "chain": chain_id,
            "res_seq": res_id[1],
            "icode": res_id[2].strip(),
            "resname": res.get_resname().strip().upper(),
            "min_dist_A": round(d, 3),
        })

    scored.sort(key=lambda x: x["min_dist_A"])

    # Apply cutoff; relax if necessary
    for limit in [cutoff, cutoff * 1.5, 30.0]:
        candidates = [s for s in scored if s["min_dist_A"] <= limit]
        if len(candidates) >= 2:
            return candidates[:2]

    # Last resort: return the two closest regardless of distance
    return scored[:2]


# ---------------------------------------------------------------------------
# Reference-structure calibration (optional)
# ---------------------------------------------------------------------------

def calibrate_from_reference(ref_pdb, carb_names):
    """
    Parse the reference PDB, find carbohydrate atoms and the two
    nearest acidic residues, and report their distances.
    Useful to choose an appropriate cutoff for the batch.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ref", ref_pdb)
    carb_atoms = collect_carb_atoms(structure, carb_names)
    if not carb_atoms:
        print(f"[WARNING] No carbohydrate atoms found in reference {ref_pdb}.")
        return None

    print(f"[Reference] {len(carb_atoms)} carbohydrate heavy atoms found.")
    results = find_two_closest_acidic(structure, carb_atoms, cutoff=30.0)
    for r in results:
        print(
            f"  Chain {r['chain']}  Res {r['resname']}{r['res_seq']}"
            f"  dist={r['min_dist_A']:.2f} Å"
        )
    if results:
        return max(r["min_dist_A"] for r in results) + 3.0   # add 3 Å buffer
    return 12.0


# ---------------------------------------------------------------------------
# Main batch processor
# ---------------------------------------------------------------------------


def find_two_closest_acidic_raw(structure, carb_coords, cutoff=12.0,
                                max_pair_dist=12.0):
    """
    Find the two acidic residues closest to the carbohydrate, subject to a
    maximum inter-residue distance constraint.

    Each residue is scored by its minimum Euclidean distance from any of its
    functional side-chain oxygens to any carbohydrate heavy atom.  Using the
    minimum-to-any-atom distance (rather than a distance to the carb centroid)
    ensures that residues near one end of a large carbohydrate are not
    penalised, and that the fallback path always returns the two residues that
    are genuinely closest to the ligand.

    The inter-residue distance filter ensures the two selected residues are
    spatially adjacent (as expected for a catalytic dyad).

    Parameters
    ----------
    cutoff : float
        Max distance in Å from any carb atom to the acidic side-chain oxygen.
        Auto-relaxed in steps if fewer than 2 candidate pairs are found.
    max_pair_dist : float
        Max distance in Å allowed between the two selected residues.
        Catalytic dyads are typically 5-10 Å apart. Default 12 Å.
    """
    acidic = collect_acidic_residues(structure)
    if not acidic:
        return []

    carb_arr = np.array(carb_coords)   # shape (N_carb, 3)

    scored = []
    for res, chain_id, res_id in acidic:
        f_atoms = functional_atoms(res)
        # shape (N_func, 3)
        f_coords = np.array([a.get_vector().get_array() for a in f_atoms])

        # ── KEY FIX ──────────────────────────────────────────────────────────
        # Score by minimum Euclidean distance from any functional oxygen to
        # any carb atom.  Previously the code used distance-to-carb-centroid,
        # which caused two residues equidistant from the centroid but on
        # opposite sides of the ligand to be selected.
        # Broadcasting: (N_func, 1, 3) - (1, N_carb, 3) → (N_func, N_carb, 3)
        dists = np.sqrt(np.sum(np.pow(f_coords[:, None, :] - carb_arr[None, :,:], 2),axis=2))                                   # shape (N_func, N_carb)
    
        d = float(dists.min())              # closest oxygen ↔ closest carb atom

        # The functional oxygen that is closest to the carb (used for the
        # inter-residue pair-distance check below)
        best_oxy_idx = int(np.unravel_index(dists.argmin(), dists.shape)[0])
        closest_fc = f_coords[best_oxy_idx]
        # ─────────────────────────────────────────────────────────────────────

        scored.append({
            "chain": chain_id,
            "res_seq": res_id[1],
            "icode": res_id[2].strip(),
            "resname": res.get_resname().strip().upper(),
            "min_dist_A": round(d, 3),
            "_coord": closest_fc,
        })

    scored.sort(key=lambda x: x["min_dist_A"])

    def _strip(entry):
        return {k: v for k, v in entry.items() if k != "_coord"}

    def _try_pairs(candidates, pair_limit):
        """Return the first pair (ranked by proximity to carb) whose
        inter-residue distance is ≤ pair_limit."""
        for i in range(len(candidates)):
            for j in range(i + 1, len(candidates)):
                pair_d = float(np.sqrt(np.sum(np.pow(candidates[i]["_coord"] - candidates[j]["_coord"],2))))
                if pair_d <= pair_limit:
                    r1, r2 = _strip(candidates[i]), _strip(candidates[j])
                    r1["pair_dist_A"] = round(pair_d, 3)
                    r2["pair_dist_A"] = round(pair_d, 3)
                    return [r1, r2]
        return None

    # ── Pass 1: honour both the distance cutoff AND max_pair_dist ────────────
    for limit in [cutoff, cutoff * 1.5, 30.0]:
        candidates = [s for s in scored if s["min_dist_A"] <= limit]
        result = _try_pairs(candidates, max_pair_dist)
        if result:
            return result

    # ── Pass 2: cutoff already at 30 Å; relax pair constraint progressively ─
    # This keeps the two residues closest to the carb while accommodating
    # dyads that are slightly further apart than max_pair_dist.
    # We do NOT fall back to ignoring the pair constraint entirely, because
    # that is what caused the "one residue near the carb, one on the far side
    # of the protein" bug when the centroid-based scoring was used.
    candidates = scored   # all residues, sorted by proximity to carb
    
    print(scored)

    for pair_limit in [max_pair_dist * 1.5, max_pair_dist * 2.0, max_pair_dist * 4.0]:
        result = _try_pairs(candidates, pair_limit)
        if result:
            return result

    # ── Absolute last resort: two closest to carb regardless of pair dist ────
    # At this point the structure is unusual; we return the best we have and
    # flag it via the large pair_dist_A value in the output.
    result = scored[:2]
    if len(result) == 2:
        pair_d = float(np.sqrt(np.sum(np.pow(result[0]["_coord"] - result[1]["_coord"], 2))))
        for r in result:
            r["pair_dist_A"] = round(pair_d, 3)
    for r in result:
        r.pop("_coord", None)
    return result


def process_batch(pdb_dir, carb_names, cutoff, max_pair_dist, out_csv):
    pdb_files = sorted(Path(pdb_dir).glob("*.pdb"))
    if not pdb_files:
        pdb_files = sorted(Path(pdb_dir).glob("*.PDB"))
    if not pdb_files:
        sys.exit(f"No .pdb files found in {pdb_dir}")

    print(f"Found {len(pdb_files)} PDB files. Processing...")

    parser = PDBParser(QUIET=True)
    rows = []

    for pdb_path in tqdm(pdb_files, unit="file"):
        row = {"pdb_file": pdb_path.name,
               "carb_atoms_found": 0,
               "res1_chain": "", "res1_name": "", "res1_seq": "", "res1_dist": "",
               "res2_chain": "", "res2_name": "", "res2_seq": "", "res2_dist": "",
               "pair_dist_A": "",
               "status": "OK"}
        try:
            # Use raw line parsing for carb detection — Biopython misparsed
            # residue names with non-alphanumeric chars (e.g. "l:b")
            carb_coords = collect_carb_atoms_raw(str(pdb_path), carb_names)
            row["carb_atoms_found"] = len(carb_coords)

            if not carb_coords:
                row["status"] = "NO_CARB"
                rows.append(row)
                continue

            # Biopython still used for acidic residue parsing (standard names only)
            structure = parser.get_structure(pdb_path.stem, str(pdb_path))
            hits = find_two_closest_acidic_raw(structure, carb_coords, cutoff, max_pair_dist)

            if len(hits) < 2:
                row["status"] = f"ONLY_{len(hits)}_ACIDIC_FOUND"

            for i, hit in enumerate(hits[:2], start=1):
                row[f"res{i}_chain"] = hit["chain"]
                row[f"res{i}_name"]  = hit["resname"]
                row[f"res{i}_seq"]   = str(hit["res_seq"]) + hit["icode"]
                row[f"res{i}_dist"]  = hit["min_dist_A"]
            if len(hits) == 2:
                row["pair_dist_A"] = hits[0].get("pair_dist_A", "")

        except Exception as exc:
            row["status"] = f"ERROR: {exc}"

        rows.append(row)

    df = pd.DataFrame(rows, columns=[
        "pdb_file", "carb_atoms_found",
        "res1_chain", "res1_name", "res1_seq", "res1_dist",
        "res2_chain", "res2_name", "res2_seq", "res2_dist",
        "pair_dist_A",
        "status",
    ])
    df.to_csv(out_csv, index=False)
    print(f"\nDone. Results written to: {out_csv}")

    # Summary
    ok  = (df["status"] == "OK").sum()
    nc  = (df["status"] == "NO_CARB").sum()
    err = len(df) - ok - nc
    print(f"  {ok:>5} structures: both residues found")
    print(f"  {nc:>5} structures: no carbohydrate detected")
    print(f"  {err:>5} structures: other issues (see 'status' column)")
    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Find two catalytic acidic residues near a carbohydrate in PDB files."
    )
    p.add_argument("--pdb_dir",   required=True,
                   help="Directory containing .pdb files")
    p.add_argument("--out",       default="acidic_residues.csv",
                   help="Output CSV file path (default: acidic_residues.csv)")
    p.add_argument("--reference", default=None,
                   help="Reference PDB file for calibrating the search cutoff")
    p.add_argument("--carb_names", default="",
                   help="Comma-separated 3-letter carbohydrate residue names "
                        "(e.g. MAN,FUC). If omitted, uses built-in sugar list.")
    p.add_argument("--cutoff",   type=float, default=12.0,
                   help="Max distance (Å) from any carb atom to acidic side chain "
                        "(default: 12.0; auto-expanded if too few hits)")
    p.add_argument("--max_pair_dist", type=float, default=12.0,
                   help="Max distance (Å) allowed between the two selected residues "
                        "(default: 12.0; catalytic dyads are typically 5-10 Å apart)")
    return p.parse_args()


def main():
    args = parse_args()
    carb_names = {n.strip().upper() for n in args.carb_names.split(",") if n.strip()}

    cutoff = args.cutoff
    if args.reference:
        print(f"Calibrating search radius using reference: {args.reference}")
        calibrated = calibrate_from_reference(args.reference, carb_names)
        if calibrated:
            cutoff = calibrated
            print(f"Using calibrated cutoff: {cutoff:.1f} Å\n")

    process_batch(args.pdb_dir, carb_names, cutoff, args.max_pair_dist, args.out)


if __name__ == "__main__":
    main()
