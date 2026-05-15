#include "md/PDBReader.h"
#include "md/PeriodicTable.h"

#include <iostream>
#include <iomanip>
#include <numeric>
#include <string>

int main(int argc, char* argv[]) {
    // Accept a PDB path on the command line, or fall back to a default.
    const std::string pdbPath = (argc > 1) ? argv[1] : "input.pdb";

    std::vector<cactus::pdb::PDBAtom> atoms;
    try {
        atoms = cactus::pdb::parsePDB(pdbPath);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    // ── Print a summary table ─────────────────────────────────────────────────
    std::cout << std::left
              << std::setw(6)  << "Serial"
              << std::setw(6)  << "Name"
              << std::setw(6)  << "Res"
              << std::setw(6)  << "Chain"
              << std::setw(7)  << "ResSeq"
              << std::setw(10) << "X(Å)"
              << std::setw(10) << "Y(Å)"
              << std::setw(10) << "Z(Å)"
              << std::setw(10) << "Occ"
              << std::setw(10) << "Mass(u)"
              << "\n"
              << std::string(79, '-') << "\n";

    std::cout << std::fixed << std::setprecision(3);

    for (const auto& a : atoms) {
        std::cout << std::left
                  << std::setw(6)  << a.serial
                  << std::setw(6)  << a.name
                  << std::setw(6)  << a.resName
                  << std::setw(6)  << a.chainID
                  << std::setw(7)  << a.resSeq
                  << std::setw(10) << a.pos.x
                  << std::setw(10) << a.pos.y
                  << std::setw(10) << a.pos.z
                  << std::setw(10) << a.occupancy
                  << std::setw(10) << a.mass
                  << "\n";
    }

    // ── Compute total mass ────────────────────────────────────────────────────
    double totalMass = std::accumulate(atoms.begin(), atoms.end(), 0.0,
        [](double sum, const cactus::pdb::PDBAtom& a) { return sum + a.mass; });

    std::cout << "\nTotal atoms : " << atoms.size()  << "\n"
              << "Total mass  : " << std::fixed << std::setprecision(3)
              << totalMass << " u\n";

    // ── Demonstrate PeriodicTable directly ───────────────────────────────────
    std::cout << "\n-- PeriodicTable direct lookup --\n";
    for (const char* sym : {"C", "N", "O", "H", "S", "Fe", "Zn"}) {
        const auto& lut = cactus::periodic::lookup();
        auto it = lut.find(sym);
        if (it != lut.end())
            std::cout << std::setw(4) << it->second->symbol
                      << "  " << std::setw(14) << it->second->name
                      << "  Z=" << std::setw(3) << it->second->atomicNumber
                      << "  mass=" << it->second->mass << " u\n";
    }

    return 0;
}
