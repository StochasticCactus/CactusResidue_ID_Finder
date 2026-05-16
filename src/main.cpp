/**
 * main.cpp
 *
 * Scans a directory of PDB files and reports pairs of adjacent acidic residues
 * (ASP / GLU) whose catalytic oxygens lie within a configurable distance of a
 * muramic-acid (or other known) carbohydrate residue.
 *
 * Results are written both to stdout and to a log file in a grep-friendly
 * format:
 *
 *   PAIR file=<path>  res1=<seqID>  res2=<seqID>  dist=<Å>
 *
 * so a full result set can be recovered with:
 *
 *   grep '^PAIR' results.log
 */

#include "md/PDBReader.h"
#include "md/PeriodicTable.h"
#include "md/Residues.h"
#include "utility/fileHandler.h"
#include "md/sorting.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdlib>

using namespace cactus;
using namespace cactus::pdb;
using namespace cactus::residues;
using namespace cactus::sort;

// ── Runtime configuration ─────────────────────────────────────────────────────

struct Config {
    std::string inputPath  = "./";       // directory that holds the PDB files
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
        << "                          (default: GLU,ASP)\n\n"
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

// ── Argument parser ───────────────────────────────────────────────────────────

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

// ── Entry point ───────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
{
    const Config cfg = parseArgs(argc, argv);

    // Full table of supported acidic residues and their catalytic oxygens.
    // Keeping this here rather than in the header makes it easy to extend.
    const std::unordered_map<std::string, std::vector<std::string>> ALL_ACIDIC_OXY = {
        {"GLU", {"OE1", "OE2"}},
        {"ASP", {"OD1", "OD2"}}
    };

    // Build the filtered map from the user-selected residues.
    std::unordered_map<std::string, std::vector<std::string>> acidicOxy;
    for (const auto& res : cfg.acidicResidues) {
        if (ALL_ACIDIC_OXY.count(res)) {
            acidicOxy[res] = ALL_ACIDIC_OXY.at(res);
        } else {
            std::cerr << "Warning: '" << res
                      << "' is not in the ACIDIC_OXY table — skipping.\n";
        }
    }

    if (acidicOxy.empty()) {
        std::cerr << "Error: no valid acidic residues selected.  Aborting.\n";
        return 1;
    }

    // Open the log file.
    std::ofstream logFile(cfg.outputFile);
    if (!logFile) {
        std::cerr << "Error: could not open output file: " << cfg.outputFile << "\n";
        return 1;
    }

    // Write a human-readable header so the log is self-describing.
    logFile << "# Acidic-residue pairs near carbohydrate residues\n"
            << "# carb-cutoff : " << cfg.carbCutoff << " Å\n"
            << "# pair-cutoff : " << cfg.pairCutoff << " Å\n"
            << "# residues    : ";
    for (size_t k = 0; k < cfg.acidicResidues.size(); ++k) {
        logFile << cfg.acidicResidues[k];
        if (k + 1 < cfg.acidicResidues.size()) logFile << ", ";
    }
    logFile << "\n#\n";

    // ── Process each PDB file in the directory ────────────────────────────────

    StructuresDir wdir(cfg.inputPath);

    for (const auto& file : wdir.dirFiles()) {

        std::vector<PDBAtom> atoms;
        try {
            atoms = parsePDB(file.path());
        } catch (const std::exception& ex) {
            std::cerr << "Error parsing " << file.path() << ": " << ex.what() << "\n";
            continue;   // skip broken files rather than aborting the whole run
        }

        // 1. Collect all atoms that belong to known carbohydrate residues.
        auto carbAtoms = CollectCarbohydrateATOMS(atoms, KNOWN_CARBS);

        // 2. Find acidic-residue oxygens within carbCutoff of any carb atom.
        auto candidates = FilterAcidicResidues(atoms, carbAtoms, acidicOxy, cfg.carbCutoff);

        // 3. Find pairs of candidate residues within pairCutoff of each other,
        //    deduplicated so each {res1, res2} pair appears only once.
        auto pairs = FilterResiduePairs(candidates, cfg.pairCutoff);

        // ── Per-file output ───────────────────────────────────────────────────

        std::cout << "### " << file.path() << " ###\n";
        logFile   << "# FILE: " << file.path() << "\n";

        if (pairs.empty()) {
            std::cout << "  (no pairs found)\n";
            logFile   << "  (no pairs found)\n";
            continue;
        }

        for (const auto& p : pairs) {
            // Format each result as a single, grep-able PAIR line.
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(3);
            oss << "PAIR"
                << "  file=" << file.path()
                << "  res1=" << p.first.first
                << "  res2=" << p.first.second
                << "  dist=" << p.second;

            const std::string line = oss.str();
            std::cout << line << "\n";
            logFile   << line << "\n";
        }
    }

    std::cout << "\nResults written to: " << cfg.outputFile << "\n";
    return 0;
}
