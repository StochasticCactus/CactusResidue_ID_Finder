/**
 * sorting.cpp
 *
 * Distance calculation and two-stage residue-pair filtering.
 *
 * Periodic boundary conditions are intentionally NOT applied.  This tool
 * targets static crystal structures and RFDiffusion models; PBC would
 * introduce spurious short-distance artefacts in those contexts.
 */

#include "md/sorting.h"
#include "md/PDBReader.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace cactus::sort {

// ── Euclidean distance ────────────────────────────────────────────────────────

double distance(const cactus::pdb::Vec3& a, const cactus::pdb::Vec3& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// ── Stage 1: filter acidic residues by proximity to a reference atom set ──────

std::vector<FAtom> FilterAcidicResidues(
        std::vector<cactus::pdb::PDBAtom>& atoms,
        std::vector<cactus::pdb::PDBAtom>& refAtoms,
        const std::unordered_map<std::string,
              std::array<std::string, 2>>& acidicOxy,
        double cutoff)
{
    std::vector<FAtom> candidates;

    cactus::pdb::Vec3 avgCoord{};
    int count = 0;

    for (const auto& atom : atoms) {

        if (count >= 2) {
            avgCoord = {};
            count    = 0;
        }

        const auto it = acidicOxy.find(atom.resName);
        if (it == acidicOxy.end())
            continue;

        const auto& oxyNames = it->second;

        if (atom.name == oxyNames[0] || atom.name == oxyNames[1]) {
            avgCoord.x += atom.pos.x * 0.5;
            avgCoord.y += atom.pos.y * 0.5;
            avgCoord.z += atom.pos.z * 0.5;
            ++count;
        }

        if (count == 2) {
            for (const auto& ref : refAtoms) {
                const double d = distance(avgCoord, ref.pos);
                if (d < cutoff) {
                    candidates.push_back({atom, d});
                }
            }
        }
    }

    return candidates;
}

// ── Stage 2: find deduplicated pairs within cutoff ────────────────────────────

std::vector<std::pair<std::pair<int,int>, double>>
FilterResiduePairs(std::vector<FAtom> atoms, double cutoff)
{
    std::map<std::pair<int,int>, double> bestDist;

    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {

            if (atoms[i].atom.resSeq == atoms[j].atom.resSeq)
                continue;

            const double d = distance(atoms[i].atom.pos, atoms[j].atom.pos);
            if (d >= cutoff)
                continue;

            const auto key = std::make_pair(
                std::min(atoms[i].atom.resSeq, atoms[j].atom.resSeq),
                std::max(atoms[i].atom.resSeq, atoms[j].atom.resSeq));

            auto it = bestDist.find(key);
            if (it == bestDist.end() || d < it->second)
                bestDist[key] = d;
        }
    }

    std::vector<std::pair<std::pair<int,int>, double>> pairs;
    pairs.reserve(bestDist.size());
    for (const auto& kv : bestDist)
        pairs.push_back({kv.first, kv.second});

    std::sort(pairs.begin(), pairs.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });

    return pairs;
}

} // namespace cactus::sort
