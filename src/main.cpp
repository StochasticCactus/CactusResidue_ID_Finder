#include "md/PDBReader.h"
#include "md/PeriodicTable.h"
#include "utility/fileHandler.h"
#include <iostream>
#include <iomanip>
#include <numeric>
#include <string>

using namespace cactus;

int main(int argc, char* argv[]) {
    // Accept a PDB path on the command line, or fall back to a default.
  
    const std::string Path = (argc > 1) ? argv[1] : "./";
    StructuresDir wdir(Path);
    
    for (const auto& file: wdir.dirFiles())
    { 

      std::vector<cactus::pdb::PDBAtom> atoms;
      try {
          atoms = cactus::pdb::parsePDB(file.path());
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
    }

    return 0;
}
