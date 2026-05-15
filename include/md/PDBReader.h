#pragma once

/*
 * PDBReader.h
 *
 * Parses ATOM and HETATM records from a standard PDB file into a flat
 * vector of PDBAtom structs.  Atom masses are assigned automatically
 * using PeriodicTable.h (elementMassFromPDBName).
 *
 * Usage:
 *   #include "md/PDBReader.h"
 *
 *   auto atoms = cactus::pdb::parsePDB("structure.pdb");
 *   for (const auto& a : atoms)
 *       std::cout << a.name << "  " << a.mass << "\n";
 */

#include <string>
#include <vector>

namespace cactus::pdb {

struct Vec3 {
    double x, y, z;
};

struct PDBAtom {
    int         serial;    // atom serial number (1-based, as in PDB)
    std::string name;      // atom name, e.g. "CA"
    std::string resName;   // residue name, e.g. "ALA"
    char        chainID;   // chain identifier character
    int         resSeq;    // residue sequence number
    Vec3        pos;       // coordinates in Ångström (PDB native)
    double      occupancy; // occupancy factor (defaults to 1.0 if absent)
    double      mass;      // atomic mass (u / Da), filled by parsePDB()
};

/**
 * Parse a PDB file and return all ATOM / HETATM records.
 *
 * @param path  Path to the .pdb file.
 * @return      Vector of PDBAtom (never empty on success).
 * @throws std::runtime_error  If the file cannot be opened or contains
 *                             no ATOM/HETATM records.
 */
std::vector<PDBAtom> parsePDB(const std::string& path);

} // namespace cactus::pdb
