#include "md/PDBReader.h"
#include "md/PeriodicTable.h"
#include "md/Residues.h"
#include "utility/fileHandler.h"
#include "md/sorting.h"
#include <iostream>
#include <iomanip>
#include <numeric>
#include <string>

using namespace cactus;
using namespace cactus::pdb;
using namespace cactus::residues;
using namespace cactus::sort;

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

       auto CarbATOMS  = CollectCarbohydrateATOMS(atoms, KNOWN_CARBS);
       auto candidates = FilterAcidicResidues(atoms, CarbATOMS, ACIDIC_OXY, 8.0);
       auto pairs      = FilterResiduePairs(candidates, 20.0);
        
       std::cout << "### " << file.path() << "###" << std::endl;

       for (const auto & p : pairs)
       {
          std::cout << std::to_string(p.first.first) << "-" << std::to_string(p.first.second) << "--" 
         << std::to_string(p.second) << std::endl;
       }




        /*
        // ── print a summary table ─────────────────────────────────────────────────
        std::cout << std::left
                  << std::setw(6)  << "serial"
                  << std::setw(6)  << "name"
                  << std::setw(6)  << "res"
                  << std::setw(6)  << "chain"
                  << std::setw(7)  << "resseq"
                  << std::setw(10) << "x(å)"
                  << std::setw(10) << "y(å)"
                  << std::setw(10) << "z(å)"
                  << std::setw(10) << "occ"
                  << std::setw(10) << "mass(u)"
                  << "\n"
                  << std::string(79, '-') << "\n";

        std::cout << std::fixed << std::setprecision(3);
      
        for (const auto& a : CarbATOMS) {
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
        }*/
      
    }
    return 0;
}
