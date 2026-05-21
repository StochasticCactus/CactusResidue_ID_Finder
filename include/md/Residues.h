#pragma once
#ifndef RESIDUES_H
#define RESIDUES_H

#include "md/PDBReader.h"
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace cactus::residues {

// ── Carbohydrate residue names ────────────────────────────────────────────────

const std::vector<std::string> KNOWN_CARBS = {
    // Hexoses
    "GLC", "BGC", "BGLC", "GLA", "GAL", "MAN", "BMA", "FRU", "TAG", "SOR",
    // Pentoses
    "RIB", "ARA", "XYL", "LYX",
    // Deoxy sugars
    "FUC", "FCA", "RHA", "RHM",
    // GlcNAc / GalNAc / ManNAc
    "NAG", "NDG", "NGA", "A2G", "BNG", "MNA",
    // Neuraminic acids
    "SIA", "SLB", "NGC",
    // Disaccharides
    "LAT", "MAL", "CEL", "TRE",
    // Uronic acids
    "GCU", "IDR",
    // Other
    "SGN", "SHB", "GCS", "LSB", "LB"
};

// ── Acidic-residue carboxyl oxygens ──────────────────────────────────────────
//
// Maps each supported acidic residue to the two atom names that form its
// carboxylate group.  These are the atoms whose midpoint is used as the
// geometric probe throughout the entire pipeline (both 'dist' and 'coord').
//
// array<string,2> is used (not vector) because the count is always exactly
// two and this avoids a heap allocation per lookup.

const std::unordered_map<std::string, std::array<std::string, 2>> ACIDIC_OXY = {
    {"ASP", {"OD1", "OD2"}},
    {"GLU", {"OE1", "OE2"}}
};

// ── Functions ─────────────────────────────────────────────────────────────────

/// Case-insensitive check: is `residue` in `knownCarbs`?
bool is_carbohydrate(std::string residue,
                     const std::vector<std::string>& knownCarbs);

/// Returns the subset of `atoms` that belong to a recognised carbohydrate.
std::vector<cactus::pdb::PDBAtom> CollectCarbohydrateATOMS(
        std::vector<cactus::pdb::PDBAtom>& atoms,
        const std::vector<std::string>&    knownCarbs);

} // namespace cactus::residues

#endif // RESIDUES_H
