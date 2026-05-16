/**
 * Residues.h
 *
 * Constants and utilities for carbohydrate residue identification.
 *
 * The central lookup table `KNOWN_CARBS` covers the most common sugar codes
 * found in PDB structures, grouped by sugar class.  `ACIDIC_OXY` maps each
 * supported acidic residue to the names of its two catalytic oxygens — the
 * atoms used as the geometric probe in Stage 1 of the proximity search.
 */

#pragma once
#ifndef Residues
#define Residues

#include <string>
#include <vector>
#include <unordered_map>
#include "md/PDBReader.h"

namespace cactus::residues {

// ── Carbohydrate residue names ────────────────────────────────────────────────

const std::vector<std::string> KNOWN_CARBS = {
    // Hexoses
    "GLC", "BGC", "BGLC", "GLA", "GAL", "MAN", "BMA",
    "FRU", "TAG", "SOR",
    // Pentoses
    "RIB", "ARA", "XYL", "LYX",
    // Deoxy sugars
    "FUC", "FCA", "RHA", "RHM",
    // GlcNAc / GalNAc / ManNAc
    "NAG", "NDG", "NGA", "A2G", "BNG", "MNA",
    // Neuraminic acids
    "SIA", "SLB", "NGC",
    // Disaccharides (sometimes appear as single residues in PDB files)
    "LAT", "MAL", "CEL", "TRE",
    // Uronic acids
    "GCU", "IDR",
    // Generic / other
    "SGN", "SHB", "GCS", "l:b", "LSB", "LB"
};

// ── Acidic-residue catalytic oxygens ─────────────────────────────────────────

/**
 * Maps each supported acidic residue to the two atom names that form its
 * carboxylate group.  These are the atoms whose midpoint is compared to
 * carbohydrate atoms in FilterAcidicResidues().
 *
 * To add further residues (e.g. LYS, HIS), extend this map and add the
 * corresponding oxygen names.  The runtime `--residues` flag then controls
 * which entries are actually used.
 */
const std::unordered_map<std::string, std::vector<std::string>> ACIDIC_OXY = {
    {"GLU", {"OE1", "OE2"}},
    {"ASP", {"OD1", "OD2"}}
};

// ── Functions ─────────────────────────────────────────────────────────────────

/**
 * Case-insensitive check: is `residue` in `KNOWN_CARBS`?
 */
bool is_carbohydrate(std::string residue, std::vector<std::string> KNOWN_CARBS);

/**
 * Returns the subset of `atoms` that belong to a recognised carbohydrate
 * residue (i.e. whose resName is in `KNOWN_CARBS`).
 */
std::vector<cactus::pdb::PDBAtom> CollectCarbohydrateATOMS(
        std::vector<cactus::pdb::PDBAtom>& atoms,
        const std::vector<std::string> KNOWN_CARBS);

} // namespace cactus::residues

#endif // Residues
