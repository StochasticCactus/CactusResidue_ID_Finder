/**
 * cmd_dist.cpp
 *
 * Implementation of the 'dist' subcommand.
 * This is the original pipeline moved verbatim out of main(), with the
 * only structural change being that it lives in cactus::cmd::runDist().
 *
 * Config and argument parsing come from commandline.h (cactus::utils).
 */

#include "utility/cmd_dist.h"
#include "utility/commandline.h"
#include "utility/fileHandler.h"
#include "utility/progressBar.h"
#include "md/PDBReader.h"
#include "md/Residues.h"
#include "md/sorting.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cactus::pdb;
using namespace cactus::residues;
using namespace cactus::sort;
using namespace cactus::utils;

namespace cactus::cmd {

int runDist(int argc, char* argv[])
{
    const Config cfg = parseArgs(argc, argv);

    // Full table of supported acidic residues and their catalytic oxygens.
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

    std::ofstream logFile(cfg.outputFile);
    if (!logFile) {
        std::cerr << "Error: could not open output file: " << cfg.outputFile << "\n";
        return 1;
    }

    // Self-describing log header.
    logFile << "# [dist] Acidic-residue pairs near carbohydrate residues\n"
            << "# carb-cutoff : " << cfg.carbCutoff << " Å\n"
            << "# pair-cutoff : " << cfg.pairCutoff << " Å\n"
            << "# residues    : ";
    for (size_t k = 0; k < cfg.acidicResidues.size(); ++k) {
        logFile << cfg.acidicResidues[k];
        if (k + 1 < cfg.acidicResidues.size()) logFile << ", ";
    }
    logFile << "\n";
    if (!cfg.extraCarbs.empty()) {
        logFile << "# extra-carbs : ";
        for (size_t k = 0; k < cfg.extraCarbs.size(); ++k) {
            logFile << cfg.extraCarbs[k];
            if (k + 1 < cfg.extraCarbs.size()) logFile << ", ";
        }
        logFile << "\n";
    }
    logFile << "#\n";

    StructuresDir wdir(cfg.inputPath);
    const std::vector<std::filesystem::directory_entry> files(
        wdir.dirFiles(), std::filesystem::directory_iterator{});
    const size_t total = files.size();

    for (size_t idx = 0; idx < total; ++idx) {
        const auto& file = files[idx];
        cactus::util::printProgress(idx, total, file.path());

        std::vector<PDBAtom> atoms;
        try {
            atoms = parsePDB(file.path());
        } catch (const std::exception& ex) {
            std::cerr << "\nError parsing " << file.path()
                      << ": " << ex.what() << "\n";
            continue;
        }

        auto knownCarbs = KNOWN_CARBS;
        knownCarbs.insert(knownCarbs.end(), cfg.extraCarbs.begin(), cfg.extraCarbs.end());

        auto carbAtoms  = CollectCarbohydrateATOMS(atoms, knownCarbs);
        auto candidates = FilterAcidicResidues(atoms, carbAtoms, acidicOxy, cfg.carbCutoff);
        auto pairs      = FilterResiduePairs(candidates, cfg.pairCutoff);

        logFile << "# FILE: " << file.path() << "\n";

        if (pairs.empty()) {
            logFile << "  (no pairs found)\n";
            continue;
        }

        for (const auto& p : pairs) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(3);
            oss << "PAIR"
                << "  file=" << file.path()
                << "  res1=" << p.first.first
                << "  res2=" << p.first.second
                << "  dist=" << p.second;
            logFile << oss.str() << "\n";
        }
    }

    cactus::util::printProgress(total, total);
    std::cout << "Results written to: " << cfg.outputFile << "\n";
    return 0;
}

} // namespace cactus::cmd
