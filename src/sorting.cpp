#include "md/sorting.h"
#include "math.h"
#include "md/PDBReader.h"
#include <algorithm>

namespace cactus::sort {
    
  double distance(const Vec3& a, const Vec3& b)
  {  
    //Intentionally NOT using periodic boundary conditions here!!!
    //The use case for which this is meant doesn't use a simulation box!!!! 
    //PBC would generate artificats!
      return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2) + std::pow(a.z - b.z, 2));
  }

  std::vector<FAtom> FilterAcidicResidues(std::vector<PDBAtom> & atoms, std::vector<PDBAtom> & CarbAtoms, const std::unordered_map<std::string, std::vector<std::string>> ACIDIC_OXY, double cutoff)
  {
      Vec3 avgCoord;
      int count = 0;
      double dist;
      FAtom candidate;
      std::vector<FAtom> candidates;

      for (const auto & atom : atoms)
      {
        if (count >= 2)
        {  
          avgCoord.x = 0;
          avgCoord.y = 0;
          avgCoord.z = 0;
          count = 0;
        }

        if (ACIDIC_OXY.count(atom.resName) <= 0)
        {
          continue;
        }

        if (ACIDIC_OXY.at(atom.resName)[0] == atom.name || ACIDIC_OXY.at(atom.resName)[1] == atom.name)
        {
          avgCoord.x += atom.pos.x / (double) 2.0;
          avgCoord.y += atom.pos.y / (double) 2.0;
          avgCoord.z += atom.pos.z / (double) 2.0;
          count++;
        }

        if (count == 2) 
        {  
           for (const auto & catom : CarbAtoms)
           { 
              dist = distance(avgCoord, catom.pos);

              if (dist < cutoff)
              {    
                  candidate.atom = atom;
                  candidate.dist = dist;
                  candidates.push_back(candidate);
              }
           } 
        }
      }
      
      return candidates; 
  }

  std::vector<std::pair<std::pair<int, int>, double>> FilterResiduePairs(std::vector<FAtom> atoms, double cutoff)
  {    
    double dist;
    std::pair<std::pair<int, int>, double> pairIJdist;
    std::vector<std::pair<std::pair<int, int>, double>> pairs;

    for (size_t i = 0; i < atoms.size(); i++)
    {
        for (size_t j = i + 1; j < atoms.size(); j++)
        {             
            if (atoms[i].atom.resSeq != atoms[j].atom.resSeq)  
            {
                dist = distance(atoms[i].atom.pos, atoms[j].atom.pos);
                if (dist < cutoff)
                { 
                    pairIJdist = {{atoms[i].atom.resSeq,atoms[j].atom.resSeq}, dist};
                    pairs.push_back(pairIJdist);
                }
            }
                      
        }
    }
    
    std::sort(pairs.begin(), pairs.end(),
    [](const auto& a, const auto& b) {
        return a.second < b.second; // ascending
    });

    return pairs;
  }
}
