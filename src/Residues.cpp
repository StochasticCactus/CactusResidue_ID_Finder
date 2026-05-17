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

namespace cactus::residues {

bool is_carbohydrate(std::string residue,
                     const std::vector<std::string>& knownCarbs)
{
    std::transform(residue.begin(), residue.end(), residue.begin(), ::toupper);
    return std::find(knownCarbs.begin(), knownCarbs.end(), residue)
           != knownCarbs.end();
}

std::vector<cactus::pdb::PDBAtom> CollectCarbohydrateATOMS(
        std::vector<cactus::pdb::PDBAtom>& atoms,
        const std::vector<std::string>&    knownCarbs)
{
    std::vector<cactus::pdb::PDBAtom> carbAtoms;
    for (const auto& atom : atoms)
        if (is_carbohydrate(atom.resName, knownCarbs))
            carbAtoms.push_back(atom);
    return carbAtoms;
}

} // namespace cactus::residues
