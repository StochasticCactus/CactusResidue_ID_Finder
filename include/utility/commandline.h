#pragma once
#ifndef commandline
#define commandline

#include <sstream>
#include <string>
#include <vector>
namespace cactus::utils {


  // ── Runtime configuration ─────────────────────────────────────────────────────

  struct Config {
      std::string inputPath  = "./";        // directory that holds the PDB files
      std::string outputFile = "results.log";

      // Distance (Å) from an acidic-residue oxygen to any carbohydrate atom.
      // Only acidic residues whose oxygens are within this radius are considered.
      double carbCutoff = 8.0;

      // Distance (Å) between the oxygens of two candidate acidic residues.
      // Only pairs closer than this are reported.
      double pairCutoff = 20.0;

      // Acidic residue types to search for.  Must be keys in ALL_ACIDIC_OXY.
      std::vector<std::string> acidicResidues = {"GLU", "ASP"};
  };

  // ── Helper: split a comma-separated string ────────────────────────────────────

  static std::vector<std::string> splitComma(const std::string& s)
  {
      std::vector<std::string> result;
      std::stringstream ss(s);
      std::string token;
      while (std::getline(ss, token, ','))
          if (!token.empty()) result.push_back(token);
      return result;
  }

}
#endif
