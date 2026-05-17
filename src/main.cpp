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
#include "utility/progressBar.h"
#include "md/sorting.h"
#include "utility/commandline.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdlib>

using namespace cactus;
using namespace cactus::pdb;
using namespace cactus::residues;
using namespace cactus::sort;
using namespace cactus::utils;


// ── Entry point ───────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
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

    // Collect the file list upfront so we know the total count for the bar.
    // directory_iterator is an input iterator with no .size() or operator[],
    // so we materialise it into a vector using the iterator-pair constructor.
    // The default-constructed directory_iterator{} serves as the end sentinel.
    StructuresDir wdir(cfg.inputPath);
    const std::vector<std::filesystem::directory_entry> files(
        wdir.dirFiles(), std::filesystem::directory_iterator{});
    const size_t total = files.size();

    for (size_t idx = 0; idx < total; ++idx) {
        const auto& file = files[idx];

        // Draw / update the progress bar on stderr so it doesn't interfere
        // with piped stdout output.
        cactus::util::printProgress(idx, total, file.path());

        std::vector<PDBAtom> atoms;
        try {
            atoms = parsePDB(file.path());
        } catch (const std::exception& ex) {
            // Print below the bar, then let the next iteration redraw it.
            std::cerr << "\nError parsing " << file.path()
                      << ": " << ex.what() << "\n";
            continue;
        }
        auto knownCarbs = KNOWN_CARBS;   // copy, don't modify the const
        knownCarbs.insert(knownCarbs.end(), cfg.extraCarbs.begin(), cfg.extraCarbs.end());
        
        // 1. Collect all atoms that belong to known carbohydrate residues.
        auto carbAtoms = CollectCarbohydrateATOMS(atoms, knownCarbs);

        // 2. Find acidic-residue oxygens within carbCutoff of any carb atom.
        auto candidates = FilterAcidicResidues(atoms, carbAtoms, acidicOxy, cfg.carbCutoff);

        // 3. Find pairs of candidate residues within pairCutoff of each other,
        //    deduplicated so each {res1, res2} pair appears only once.
        auto pairs      = FilterResiduePairs(candidates, cfg.pairCutoff);

        // Write results to the log file only; the bar owns the terminal line.
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

    // Final call completes the bar and moves the cursor to the next line.
    cactus::util::printProgress(total, total);

    std::cout << "Results written to: " << cfg.outputFile << "\n";
    return 0;
}
