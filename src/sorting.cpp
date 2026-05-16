#include "md/sorting.h"
#include "math.h"

namespace md::sorting
{
  double distance(Vec3 & a, Vec3 b);
  {  
    //Intentionally NOT using periodic boundary conditions here!!!
    //The use case for which this is meant doesn't use a simulation box!!!! 
    //PBC would generate artificats!
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2) + std::pow(a.z - b.z, 2)) ;
  }

    std::vector<PDBAtom> sortAcidicResidues(std::vector<PDBAtoms> & atoms, const std::unordered_map<std::string, std::vector<std::string>> ACIDIC_OXY)
    {
      for (const auto & atom : atoms)
      {

        if (ACIDIC_OXY.count(atom.residueName) > 0)
        {
          

        }         
      }
    }

}
