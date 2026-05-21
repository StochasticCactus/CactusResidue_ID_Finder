#pragma once
#ifndef SORTING_H
#define SORTING_H

#include "md/PDBReader.h"
#include "md/Residues.h"
#include <unordered_map>
#include <utility>
#include <vector>

namespace cactus::sort {

// ── Data structures ───────────────────────────────────────────────────────────

/// An acidic residue that passed the carbohydrate-proximity filter, together
/// with the distance from its carboxyl midpoint to the nearest carbohydrate atom.
struct FAtom {
    cactus::pdb::PDBAtom atom;
    double               dist;
};

// ── Functions ─────────────────────────────────────────────────────────────────

/// Euclidean distance between two 3-D points.
/// Periodic boundary conditions are intentionally not applied.
double distance(const cactus::pdb::Vec3& a, const cactus::pdb::Vec3& b);

/**
 * Stage 1 — filter acidic residues by proximity to a reference atom set.
 *
 * For each acidic residue in `atoms` (types determined by `acidicOxy`),
 * computes the midpoint of its two carboxyl oxygens and checks whether that
 * midpoint lies within `cutoff` Å of any atom in `refAtoms`.  Passing
 * residues are returned as FAtom entries.
 *
 * `refAtoms` is generic: it can be carbohydrate atoms ('dist' subcommand) or
 * any other reference set ('coord' subcommand).
 *
 * @param atoms      Full atom list from parsePDB().
 * @param refAtoms   Reference atom set to measure distance against.
 * @param acidicOxy  Map: residue name → {oxygen1_name, oxygen2_name}.
 * @param cutoff     Distance threshold in Å.
 */
std::vector<FAtom> FilterAcidicResidues(
        std::vector<cactus::pdb::PDBAtom>& atoms,
        std::vector<cactus::pdb::PDBAtom>& refAtoms,
        const std::unordered_map<std::string,
              std::array<std::string, 2>>& acidicOxy,
        double cutoff);

/**
 * Stage 2 — find deduplicated pairs of candidate residues within cutoff.
 *
 * Each pair is keyed on {min(resSeq), max(resSeq)} so the same two residues
 * are never reported more than once.  Only the minimum distance observed for
 * each unique pair is retained.
 *
 * @param atoms   Output of FilterAcidicResidues().
 * @param cutoff  Maximum inter-residue distance (Å) to report.
 * @returns       Unique pairs: {{resSeq1, resSeq2}, dist}, sorted ascending.
 */
std::vector<std::pair<std::pair<int,int>, double>>
FilterResiduePairs(std::vector<FAtom> atoms, double cutoff);

} // namespace cactus::sort

#endif // SORTING_H
