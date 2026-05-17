#pragma once
#ifndef commandline
#define commandline
#include <iostream>
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
      //Carbohydrate residue names to search for.
      std::vector<std::string> extraCarbs = {};
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

  // ── Help text ─────────────────────────────────────────────────────────────────

  static void printHelp(const char* prog)
  {
      std::cout
          << "Usage: " << prog << " [OPTIONS] [PATH]\n\n"
          << "Find pairs of adjacent acidic residues (ASP / GLU) in the vicinity\n"
          << "of muramic-acid carbohydrate residues across all PDB files in PATH.\n\n"
          << "Arguments:\n"
          << "  PATH                    Directory containing PDB files (default: ./)\n\n"
          << "Options:\n"
          << "  -h, --help              Show this help message and exit\n"
          << "  -o, --output FILE       Write results to FILE (default: results.log)\n"
          << "  --carb-cutoff DIST      Oxygen-to-carbohydrate distance cutoff in Å\n"
          << "                          (default: 8.0)\n"
          << "  --pair-cutoff DIST      Acidic-residue pair distance cutoff in Å\n"
          << "                          (default: 20.0)\n"
          << "  --residues RES          Comma-separated acidic residue types to consider\n"
          << "                          (default: GLU,ASP)\n"
          << "  --carbs RES             Comma-separated extra carbohydrate residue names\n"
          << "                          \n\n"
          << "Output format (log file and stdout):\n"
          << "  Each result is a single line beginning with 'PAIR':\n"
          << "    PAIR file=<path>  res1=<seqID>  res2=<seqID>  dist=<Å>\n\n"
          << "  Recover all pairs from the log with:\n"
          << "    grep '^PAIR' results.log\n\n"
          << "  Or filter by file:\n"
          << "    grep '^PAIR.*1abc' results.log\n\n"
          << "Examples:\n"
          << "  " << prog << " ./structures/ --carb-cutoff 6.0 --pair-cutoff 15.0\n"
          << "  " << prog << " ./structures/ -o my_run.log --residues GLU\n";
  }

  static Config parseArgs(int argc, char* argv[])
  {
      Config cfg;
   
      for (int i = 1; i < argc; ++i) {
          std::string arg = argv[i];
   
          if (arg == "-h" || arg == "--help") {
              printHelp(argv[0]);
              std::exit(0);
   
          } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
              cfg.outputFile = argv[++i];
   
          } else if (arg == "--carb-cutoff" && i + 1 < argc) {
              cfg.carbCutoff = std::stod(argv[++i]);
   
          } else if (arg == "--pair-cutoff" && i + 1 < argc) {
              cfg.pairCutoff = std::stod(argv[++i]);
   
          } else if (arg == "--residues" && i + 1 < argc) {
              cfg.acidicResidues = splitComma(argv[++i]);
   
          } else if (arg == "--carbs" && i + 1 < argc) {
              cfg.extraCarbs = splitComma(argv[++i]);
   
          } else if (arg[0] != '-') {
              cfg.inputPath = arg;
   
          } else {
              std::cerr << "Unknown option: " << arg
                        << "  (run with --help for usage)\n";
              std::exit(1);
          }
      }
      return cfg;
  }
}
#endif
