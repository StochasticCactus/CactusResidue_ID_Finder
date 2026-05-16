#include "md/Residues.h"
#include <algorithm>
#include <string>
#include "md/PDBReader.h"

using namespace cactus::pdb;

namespace cactus::residues {
  
    bool is_carbohydrate(std::string residue, const std::vector<std::string> KNOWN_CARBS)
    {
      std::transform(residue.begin(), residue.end(), residue.begin(), ::toupper);
      return std::find(KNOWN_CARBS.begin(), KNOWN_CARBS.end(), residue) != KNOWN_CARBS.end();
    }
    
    std::vector<PDBAtom> CollectCarbohydrateATOMS(std::vector<PDBAtom> & atoms, const std::vector<std::string> KNOWN_CARBS)
    {
      std::vector<PDBAtom> CarbAtoms;

      for (const auto & atom : atoms)
      {
        if (is_carbohydrate(atom.resName, KNOWN_CARBS))
          
          CarbAtoms.push_back(atom);
      }
      return CarbAtoms;
    } 

    
    std::vector<PDBAtom> sortAcidicResidues(std::vector<PDBAtoms> & atoms, const std::unordered_map<std::string, std::vector<std::string>> ACIDIC_OXY)
    {
      for (const auto & atom : atoms)
      {
        


        
      }


    }


  
}
