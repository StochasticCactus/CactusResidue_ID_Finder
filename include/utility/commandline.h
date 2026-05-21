#pragma once
#ifndef COMMANDLINE_H
#define COMMANDLINE_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace cactus::utils {

// ── Runtime configuration for the 'dist' subcommand ──────────────────────────

struct Config {
    std::string inputPath  = "./";
    std::string outputFile = "results.log";

    // Distance (Å) from an acidic-residue carboxyl midpoint to any carbohydrate
    // atom.  Only residues within this radius are considered candidates.
    double carbCutoff = 8.0;

    // Distance (Å) between two candidate carboxyl midpoints.
    // Only pairs closer than this are reported.
    double pairCutoff = 20.0;

    // Extra carbohydrate residue names supplied at runtime via --carbs.
    std::vector<std::string> extraCarbs = {};

    // Acidic residue types to search for.  Must be keys in ACIDIC_OXY.
    std::vector<std::string> acidicResidues = {"GLU", "ASP"};
};

// ── String utilities (used by both subcommands) ───────────────────────────────

/// Split a comma-separated string into tokens, discarding empty fields.
inline std::vector<std::string> splitComma(const std::string& s)
{
    std::vector<std::string> result;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ','))
        if (!token.empty()) result.push_back(token);
    return result;
}

/// Split a comma-separated string of integers.
inline std::vector<int> splitInt(const std::string& s)
{
    std::vector<int> result;
    for (const auto& tok : splitComma(s))
        result.push_back(std::stoi(tok));
    return result;
}

// ── Help text and argument parsing for 'dist' ─────────────────────────────────

inline void printDistHelp(const char* prog)
{
    std::cout
        << "Usage: " << prog << " dist [OPTIONS] [PATH]\n\n"
        << "Find pairs of adjacent acidic residues (ASP / GLU) in the vicinity\n"
        << "of carbohydrate residues across all PDB files in PATH.\n\n"
        << "Arguments:\n"
        << "  PATH                    Directory containing PDB files (default: ./)\n\n"
        << "Options:\n"
        << "  -h, --help              Show this help message and exit\n"
        << "  -o, --output FILE       Write results to FILE (default: results.log)\n"
        << "  --carb-cutoff DIST      Carboxyl-midpoint to carbohydrate distance\n"
        << "                          cutoff in Å (default: 8.0)\n"
        << "  --pair-cutoff DIST      Acidic-residue pair distance cutoff in Å\n"
        << "                          (default: 20.0)\n"
        << "  --residues RES          Comma-separated acidic residue types\n"
        << "                          (default: GLU,ASP)\n"
        << "  --carbs RES             Comma-separated extra carbohydrate residue names\n\n"
        << "Output format:\n"
        << "  PAIR  file=<path>  res1=<seqID>  res2=<seqID>  dist=<Å>\n\n"
        << "  Recover all pairs with:\n"
        << "    grep '^PAIR' results.log\n\n"
        << "Examples:\n"
        << "  " << prog << " dist ./structures/ --carb-cutoff 6.0 --pair-cutoff 15.0\n"
        << "  " << prog << " dist ./structures/ -o my_run.log --residues GLU\n";
}

inline Config parseDistArgs(int argc, char* argv[])
{
    Config cfg;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printDistHelp(argv[0]);
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
                      << "  (run '" << argv[0] << " dist --help' for usage)\n";
            std::exit(1);
        }
    }
    return cfg;
}

} // namespace cactus::utils

#endif // COMMANDLINE_H
