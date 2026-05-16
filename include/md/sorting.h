/**
 * sorting.h
 *
 * Distance calculation and two-stage residue-pair filtering for the
 * acidic-residue / carbohydrate proximity search.
 *
 * Pipeline overview:
 *
 *   FilterAcidicResidues()  →  FilterResiduePairs()
 *
 *   Stage 1 selects acidic residues (e.g. GLU / ASP) whose catalytic oxygens
 *   lie within a user-defined cutoff of any carbohydrate atom.
 *
 *   Stage 2 takes those candidates and returns deduplicated pairs of distinct
 *   residues that are within a second cutoff of each other, sorted by
 *   ascending distance.
 */

#pragma once
#ifndef sorting
#define sorting

#include "md/PDBReader.h"
#include <unordered_map>
#include <vector>
#include <utility>

using namespace cactus::pdb;

namespace cactus::sort {

// ── Data structures ───────────────────────────────────────────────────────────

/**
 * An atom that has passed the carbohydrate-proximity filter, together with
 * the distance from its midpoint coordinate to the nearest carbohydrate atom.
 */
struct FAtom {
    PDBAtom atom;   ///< The underlying PDB atom record.
    double  dist;   ///< Distance (Å) to the nearest carbohydrate atom.
};

/**
 * A pair of PDB atoms (unused by the current pipeline but kept for potential
 * downstream use, e.g. reporting the specific interacting atoms).
 */
struct AtomPair {
    PDBAtom atom1;
    PDBAtom atom2;
};

// ── Functions ─────────────────────────────────────────────────────────────────

/**
 * Euclidean distance between two 3-D points.
 *
 * Note: periodic boundary conditions are intentionally NOT applied.
 * See sorting.cpp for the rationale.
 */
double distance(const Vec3& a, const Vec3& b);

/**
 * Stage 1 — filter acidic residues by proximity to carbohydrate atoms.
 *
 * Returns one FAtom entry for each acidic-residue oxygen pair (e.g. OE1/OE2
 * of a GLU) whose midpoint lies within `cutoff` Å of any atom in `CarbAtoms`.
 *
 * @param atoms       Full atom list from parsePDB().
 * @param CarbAtoms   Carbohydrate atoms (from CollectCarbohydrateATOMS()).
 * @param ACIDIC_OXY  Map: residue name → {oxygen1_name, oxygen2_name}.
 * @param cutoff      Distance threshold in Å.
 */
std::vector<FAtom> FilterAcidicResidues(
        std::vector<PDBAtom>& atoms,
        std::vector<PDBAtom>& CarbAtoms,
        const std::unordered_map<std::string, std::vector<std::string>> ACIDIC_OXY,
        double cutoff);

/**
 * Stage 2 — find deduplicated pairs of candidate residues within cutoff.
 *
 * Each pair is identified by {min(resSeq), max(resSeq)} so that the same two
 * residues are never reported more than once, regardless of how many FAtom
 * entries were produced for each of them in Stage 1.  Only the minimum
 * distance observed for each unique pair is retained.
 *
 * @param atoms   Output of FilterAcidicResidues().
 * @param cutoff  Maximum inter-residue distance (Å) to report.
 * @returns       Unique pairs: {{resSeq1, resSeq2}, dist}, sorted ascending.
 */
std::vector<std::pair<std::pair<int,int>, double>>
FilterResiduePairs(std::vector<FAtom> atoms, double cutoff);

} // namespace cactus::sort

#endif // sorting
