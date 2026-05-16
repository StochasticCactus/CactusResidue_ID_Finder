/**
 * Residues.cpp
 *
 * Utilities for identifying and collecting carbohydrate atoms from a parsed
 * PDB atom list.
 */

#include "md/Residues.h"
#include "md/PDBReader.h"
#include <algorithm>
#include <string>

using namespace cactus::pdb;

namespace cactus::residues {

/**
 * Returns true when `residue` (case-insensitive) is present in `KNOWN_CARBS`.
 */
bool is_carbohydrate(std::string residue, const std::vector<std::string> KNOWN_CARBS)
{
    std::transform(residue.begin(), residue.end(), residue.begin(), ::toupper);
    return std::find(KNOWN_CARBS.begin(), KNOWN_CARBS.end(), residue)
           != KNOWN_CARBS.end();
}

/**
 * Extracts every atom whose residue name is listed in `KNOWN_CARBS`.
 *
 * These atoms feed Stage 1 of the pipeline (FilterAcidicResidues) as the
 * reference set against which acidic-residue oxygens are compared.
 *
 * @param atoms       Full atom list from parsePDB().
 * @param KNOWN_CARBS List of recognised carbohydrate residue names.
 * @returns           Subset of `atoms` belonging to carbohydrate residues.
 */
std::vector<PDBAtom> CollectCarbohydrateATOMS(
        std::vector<PDBAtom>& atoms,
        const std::vector<std::string> KNOWN_CARBS)
{
    std::vector<PDBAtom> carbAtoms;
    for (const auto& atom : atoms)
        if (is_carbohydrate(atom.resName, KNOWN_CARBS))
            carbAtoms.push_back(atom);
    return carbAtoms;
}

} // namespace cactus::residues
