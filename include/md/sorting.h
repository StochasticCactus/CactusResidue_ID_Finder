#pragma once
#ifndef sorting
#define sorting
#include "md/PDBReader.h"
#include <unordered_map>
using namespace cactus::pdb;

namespace cactus::sort {
  
  struct FAtom {
    
    PDBAtom atom;
    double  dist;

  };

  struct AtomPair 
  {
    PDBAtom atom1;
    PDBAtom atom2;
  };
  
  double distance(const Vec3& a, const Vec3& b);
  //Intentionally NOT using periodic boundary conditions here!!!
  
  std::vector<FAtom> FilterAcidicResidues(std::vector<PDBAtom> & atoms, std::vector<PDBAtom> & CarbAtoms, const std::unordered_map<std::string, std::vector<std::string>> ACIDIC_OXY, double cutoff);

  std::vector<std::pair<std::pair<int, int>, double>> FilterResiduePairs(std::vector<FAtom> atoms, double cutoff);
}
#endif
