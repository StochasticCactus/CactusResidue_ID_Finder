/**
 * sorting.cpp
 *
 * Distance calculations and residue-pair filtering for the acidic-residue /
 * carbohydrate proximity search.
 *
 * Two-stage pipeline:
 *
 *  1. FilterAcidicResidues  – keeps acidic residues whose catalytic oxygens
 *                             are within `cutoff` Å of any carbohydrate atom.
 *
 *  2. FilterResiduePairs    – from the surviving candidates, reports every
 *                             pair of distinct residues within a second cutoff,
 *                             deduplicated and sorted by ascending distance.
 */

#include "md/sorting.h"
#include "math.h"
#include "md/PDBReader.h"
#include <algorithm>
#include <map>

namespace cactus::sort {

// ── Euclidean distance (no PBC) ───────────────────────────────────────────────

double distance(const Vec3& a, const Vec3& b)
{
    // Periodic boundary conditions are intentionally NOT applied here.
    // The structures this tool targets (e.g. cell-wall glycopeptides in
    // crystal structures) are analysed as static, non-periodic assemblies.
    // Applying PBC would introduce spurious short-distance artefacts.
    return std::sqrt(  std::pow(a.x - b.x, 2)
                     + std::pow(a.y - b.y, 2)
                     + std::pow(a.z - b.z, 2));
}

// ── Stage 1: filter acidic residues by proximity to carbohydrate atoms ────────

/**
 * Iterates through `atoms` and identifies acidic residues (GLU / ASP, or
 * whatever is in `ACIDIC_OXY`) whose two catalytic oxygens lie within
 * `cutoff` Å of at least one atom in `CarbAtoms`.
 *
 * The function works by accumulating the positions of the two oxygens of each
 * residue (e.g. OE1 + OE2 for GLU) into a running average, then — once both
 * oxygens have been seen — comparing that average position against every
 * carbohydrate atom.  It therefore relies on OE1 and OE2 (and their ASP
 * counterparts) appearing consecutively in the atom list, which is standard
 * PDB column ordering.
 *
 * @param atoms       Full atom list from the PDB file.
 * @param CarbAtoms   Pre-filtered list of carbohydrate atoms.
 * @param ACIDIC_OXY  Map of residue name → {oxygen1, oxygen2} to match.
 * @param cutoff      Distance threshold in Å.
 * @returns           Vector of FAtom (atom + nearest-carb distance) for every
 *                    acidic-residue oxygen that passed the cutoff check.
 */
std::vector<FAtom> FilterAcidicResidues(
        std::vector<PDBAtom>& atoms,
        std::vector<PDBAtom>& CarbAtoms,
        const std::unordered_map<std::string, std::vector<std::string>> ACIDIC_OXY,
        double cutoff)
{
    Vec3 avgCoord{};    // accumulates the midpoint of the two catalytic oxygens
    int  count = 0;     // how many matching oxygens we have collected so far
    std::vector<FAtom> candidates;

    for (const auto& atom : atoms) {

        // Once we have seen both oxygens of a residue (count == 2), reset
        // ready for the next residue on the next iteration.
        if (count >= 2) {
            avgCoord = {};
            count    = 0;
        }

        // Skip atoms whose residue type is not in the acidic-oxygen table.
        if (ACIDIC_OXY.count(atom.resName) == 0)
            continue;

        const auto& oxyNames = ACIDIC_OXY.at(atom.resName);

        // Accumulate position if this is one of the two catalytic oxygens.
        if (atom.name == oxyNames[0] || atom.name == oxyNames[1]) {
            avgCoord.x += atom.pos.x / 2.0;
            avgCoord.y += atom.pos.y / 2.0;
            avgCoord.z += atom.pos.z / 2.0;
            ++count;
        }

        // Once both oxygens are in, measure to every carbohydrate atom.
        if (count == 2) {
            for (const auto& catom : CarbAtoms) {
                double dist = distance(avgCoord, catom.pos);
                if (dist < cutoff) {
                    candidates.push_back({atom, dist});
                }
            }
        }
    }

    return candidates;
}

// ── Stage 2: find pairs of candidate residues within cutoff ───────────────────

/**
 * From the candidate list produced by FilterAcidicResidues, finds every pair
 * of atoms belonging to *different* residues that are within `cutoff` Å of
 * each other.
 *
 * Deduplication:
 *   Because FilterAcidicResidues can emit several FAtom entries for the same
 *   residue (e.g. one for each nearby carbohydrate atom), the naive O(n²) loop
 *   would emit the same {res1, res2} pair multiple times.  We deduplicate by
 *   keying on the ordered pair {min(resSeq), max(resSeq)} and retaining only
 *   the minimum distance seen for that pair.
 *
 * @param atoms   Candidate acidic-residue atoms.
 * @param cutoff  Maximum inter-residue distance (Å) to report.
 * @returns       Unique pairs sorted by ascending distance.
 */
std::vector<std::pair<std::pair<int,int>, double>>
FilterResiduePairs(std::vector<FAtom> atoms, double cutoff)
{
    // std::map keeps keys ordered and makes duplicate detection trivial.
    // Key:   {min(resSeq), max(resSeq)}  — canonical ordering so that
    //         {50, 60} and {60, 50} map to the same entry.
    // Value: best (minimum) distance seen for that pair.
    std::map<std::pair<int,int>, double> bestDist;

    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {

            // Atoms from the same residue are not a "pair".
            if (atoms[i].atom.resSeq == atoms[j].atom.resSeq)
                continue;

            const double dist = distance(atoms[i].atom.pos, atoms[j].atom.pos);
            if (dist >= cutoff)
                continue;

            const auto key = std::make_pair(
                std::min(atoms[i].atom.resSeq, atoms[j].atom.resSeq),
                std::max(atoms[i].atom.resSeq, atoms[j].atom.resSeq));

            auto it = bestDist.find(key);
            if (it == bestDist.end() || dist < it->second)
                bestDist[key] = dist;
        }
    }

    // Flatten the map into the return vector and sort by distance.
    std::vector<std::pair<std::pair<int,int>, double>> pairs;
    pairs.reserve(bestDist.size());
    for (const auto& kv : bestDist)
        pairs.push_back({kv.first, kv.second});

    std::sort(pairs.begin(), pairs.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });

    return pairs;
}

} // namespace cactus::sort
