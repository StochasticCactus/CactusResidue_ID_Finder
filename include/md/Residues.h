#pragma once
#ifndef Residues
#define Residues
#include <string>
#include <vector>
#include <unordered_map>

namespace cactus::residues 
{

  std::vector<std::string> KNOWN_CARBS = {
    //Hexoses
    "GLC", "BGC", "BGLC", "GLA", "GAL", "MAN", "BMA",
    "FRU", "TAG", "SOR",
    //Pentoses
    "RIB", "ARA", "XYL", "LYX",
    //Deoxy sugars
    "FUC", "FCA", "RHA", "RHM",
    //GlcNAc / GalNAc / ManNAc
    "NAG", "NDG", "NGA", "A2G", "BNG", "MNA",
    //Neuraminic acids
    "SIA", "SLB", "NGC",
    //Disaccharides (sometimes appear as single residues)
    "LAT", "MAL", "CEL", "TRE",
    //Uronic acids
    "GCU", "IDR",
    //Generic / other
    "SGN", "SHB", "GCS", "l:b", "LSB", "LB"};

    std::unordered_map<std::string, std::vector<std::string>>  ACIDIC_OXY;
    ACIDIC_OXY["GLU"] = {"OE1", "OE2"};
    ACIDIC_OXY["ASP"] = {"OD1", "OD2"};

}

#endif
