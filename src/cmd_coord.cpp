/**
 * cmd_coord.cpp
 *
 * Implementation of the 'coord' subcommand.
 *
 * A reference PDB file is parsed first and a subset of its atoms is extracted
 * as anchor points (the "reference set").  For each query PDB in the input
 * directory the same two-stage pipeline used by 'dist' is then applied, with
 * the reference atoms playing the role that carbohydrate atoms play in 'dist':
 *
 *   Stage 1 — FilterAcidicResidues:  keep acidic residues whose carboxyl
 *             midpoint (average of the two oxygens) is within --carb-cutoff
 *             of any reference atom.
 *
 *   Stage 2 — FilterResiduePairs:    from the survivors, report deduplicated
 *             pairs within --pair-cutoff, sorted by ascending distance.
 *
 * This design means 'coord' can be used with any reference geometry — a
 * ligand, a sugar, another chain, a set of water molecules, etc. — simply by
 * choosing appropriate --ref-residues and --ref-atoms filters.
 */

#include "utility/cmd_coord.h"
#include "utility/fileHandler.h"
#include "utility/progressBar.h"
#include "md/PDBReader.h"
#include "md/sorting.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cactus::pdb;
using namespace cactus::sort;

namespace cactus::cmd {

// ── Argument parsing (private to this translation unit) ───────────────────────

static void printCoordHelp(const char* prog)
{
    std::cout
        << "Usage: " << prog << " --ref FILE [OPTIONS] [PATH]\n\n"
        << "Find pairs of acidic residues (ASP / GLU) whose carboxyl groups\n"
        << "lie within a cutoff distance of atoms in a reference PDB structure,\n"
        << "across all query PDB files in PATH.\n\n"
        << "Arguments:\n"
        << "  PATH                    Directory of query PDB files (default: ./)\n\n"
        << "Required:\n"
        << "  --ref FILE              Reference PDB file that defines the anchor\n"
        << "                          coordinate region\n\n"
        << "Options:\n"
        << "  -h, --help              Show this help message and exit\n"
        << "  -o, --output FILE       Write results to FILE (default: coord_results.log)\n"
        << "  --ref-residues RES      Comma-separated residue names to extract from the\n"
        << "                          reference (default: all residues)\n"
        << "  --ref-atoms ATOMS       Comma-separated atom names to extract from those\n"
        << "                          residues (default: all atoms)\n"
        << "  --ref-seqid IDS         Comma-separated residue sequence numbers to extract\n"
        << "                          from the reference (e.g. 142,158; default: all)\n"
        << "  --carb-cutoff DIST      Carboxyl-midpoint to reference-atom cutoff in Å\n"
        << "                          (default: 8.0)\n"
        << "  --pair-cutoff DIST      Acidic-residue pair distance cutoff in Å\n"
        << "                          (default: 20.0)\n"
        << "  --residues RES          Comma-separated acidic residue types to search for\n"
        << "                          in the query structures (default: GLU,ASP)\n\n"
        << "Output format:\n"
        << "  PAIR  file=<path>  res1=<seqID>  res2=<seqID>  dist=<Å>\n\n"
        << "  Recover all pairs with:\n"
        << "    grep '^PAIR' coord_results.log\n\n"
        << "Examples:\n"
        << "  # Use all atoms of a ligand residue LIG as the reference:\n"
        << "  " << prog << " --ref ref.pdb --ref-residues LIG ./structures/\n\n"
        << "  # Use only the carboxyl oxygens of a known ASP in the reference:\n"
        << "  " << prog << " --ref ref.pdb --ref-residues ASP --ref-atoms OD1,OD2 ./structures/\n\n"
        << "  # Wider search radius, GLU only:\n"
        << "  " << prog << " --ref ref.pdb --carb-cutoff 12.0 --residues GLU ./structures/\n";
}

static std::vector<std::string> splitComma(const std::string& s)
{
    std::vector<std::string> result;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ','))
        if (!token.empty()) result.push_back(token);
    return result;
}

static std::vector<int> splitInt(const std::string& s)
{
    std::vector<int> result;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ','))
        if (!token.empty()) result.push_back(std::stoi(token));
    return result;
}

static CoordConfig parseCoordArgs(int argc, char* argv[])
{
    CoordConfig cfg;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printCoordHelp(argv[0]);
            std::exit(0);

        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            cfg.outputFile = argv[++i];

        } else if (arg == "--ref" && i + 1 < argc) {
            cfg.refPath = argv[++i];

        } else if (arg == "--ref-residues" && i + 1 < argc) {
            cfg.refResidues = splitComma(argv[++i]);

        } else if (arg == "--ref-atoms" && i + 1 < argc) {
            cfg.refAtoms = splitComma(argv[++i]);

        } else if (arg == "--ref-seqid" && i + 1 < argc) {
            cfg.refSeqIDs = splitInt(argv[++i]);

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

// ── Reference atom extraction ─────────────────────────────────────────────────

/**
 * Parse the reference PDB and return the subset of atoms selected by the
 * --ref-residues and --ref-atoms filters.
 *
 * An empty filter list means "accept everything" for that dimension, so:
 *   refResidues = {},  refAtoms = {}   →  every atom in the reference
 *   refResidues = {"LIG"}, refAtoms = {}   →  all atoms of residue LIG
 *   refResidues = {"ASP"}, refAtoms = {"OD1","OD2"}  →  ASP carboxyl oxygens only
 */
static std::vector<PDBAtom> extractReferenceAtoms(
        const std::string& refPath,
        const std::vector<std::string>& refResidues,
        const std::vector<std::string>& refAtoms,
        const std::vector<int>& refSeqIDs)
{
    const std::vector<PDBAtom> all = parsePDB(refPath);
    std::vector<PDBAtom> selected;

    for (const auto& atom : all) {

        // Residue filter.
        if (!refResidues.empty()) {
            bool found = false;
            for (const auto& r : refResidues)
                if (atom.resName == r) { found = true; break; }
            if (!found) continue;
        }

        // Atom-name filter.
        if (!refAtoms.empty()) {
            bool found = false;
            for (const auto& a : refAtoms)
                if (atom.name == a) { found = true; break; }
            if (!found) continue;
        }

        // Sequence-number filter.
        if (!refSeqIDs.empty()) {
            bool found = false;
            for (const int id : refSeqIDs)
                if (atom.resSeq == id) { found = true; break; }
            if (!found) continue;
        }

        selected.push_back(atom);
    }

    return selected;
}

// ── Entry point ───────────────────────────────────────────────────────────────

int runCoord(int argc, char* argv[])
{
    const CoordConfig cfg = parseCoordArgs(argc, argv);

    if (cfg.refPath.empty()) {
        std::cerr << "Error: --ref <reference.pdb> is required.\n"
                  << "Run '" << argv[0] << " --help' for usage.\n";
        return 1;
    }

    // Parse the reference PDB and extract anchor atoms.
    std::vector<PDBAtom> refAtoms;
    try {
        refAtoms = extractReferenceAtoms(cfg.refPath, cfg.refResidues, cfg.refAtoms, cfg.refSeqIDs);
    } catch (const std::exception& ex) {
        std::cerr << "Error reading reference file '" << cfg.refPath
                  << "': " << ex.what() << "\n";
        return 1;
    }

    if (refAtoms.empty()) {
        std::cerr << "Error: no atoms selected from the reference file after filtering.\n"
                  << "Check --ref-residues and --ref-atoms.\n";
        return 1;
    }

    std::cout << "Reference: " << refAtoms.size()
              << " anchor atom(s) loaded from " << cfg.refPath << "\n";

    // Build the acidic-oxygen map filtered to the user-selected residue types.
    const std::unordered_map<std::string, std::vector<std::string>> ALL_ACIDIC_OXY = {
        {"GLU", {"OE1", "OE2"}},
        {"ASP", {"OD1", "OD2"}}
    };

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

    // Open log file.
    std::ofstream logFile(cfg.outputFile);
    if (!logFile) {
        std::cerr << "Error: could not open output file: " << cfg.outputFile << "\n";
        return 1;
    }

    // Self-describing log header.
    logFile << "# [coord] Acidic-residue pairs near a reference structure\n"
            << "# reference   : " << cfg.refPath << "\n"
            << "# ref-residues: "
            << (cfg.refResidues.empty() ? "(all)" : "");
    for (size_t k = 0; k < cfg.refResidues.size(); ++k) {
        logFile << cfg.refResidues[k];
        if (k + 1 < cfg.refResidues.size()) logFile << ", ";
    }
    logFile << "\n# ref-atoms   : "
            << (cfg.refAtoms.empty() ? "(all)" : "");
    for (size_t k = 0; k < cfg.refAtoms.size(); ++k) {
        logFile << cfg.refAtoms[k];
        if (k + 1 < cfg.refAtoms.size()) logFile << ", ";
    }
    if (!cfg.refSeqIDs.empty()) {
        logFile << "\n# ref-seqids  : ";
        for (size_t k = 0; k < cfg.refSeqIDs.size(); ++k) {
            logFile << cfg.refSeqIDs[k];
            if (k + 1 < cfg.refSeqIDs.size()) logFile << ", ";
        }
    }
    logFile << "\n# anchor atoms: " << refAtoms.size() << "\n"
            << "# carb-cutoff : " << cfg.carbCutoff << " Å\n"
            << "# pair-cutoff : " << cfg.pairCutoff << " Å\n"
            << "# residues    : ";
    for (size_t k = 0; k < cfg.acidicResidues.size(); ++k) {
        logFile << cfg.acidicResidues[k];
        if (k + 1 < cfg.acidicResidues.size()) logFile << ", ";
    }
    logFile << "\n#\n";

    // ── Query loop ────────────────────────────────────────────────────────────

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

        // Stage 1: find acidic residues whose carboxyl midpoint is within
        // carbCutoff of any reference atom.  We reuse FilterAcidicResidues
        // directly — refAtoms plays the same role as CarbAtoms in 'dist'.
        auto candidates = FilterAcidicResidues(atoms, refAtoms, acidicOxy, cfg.carbCutoff);

        // Stage 2: deduplicated pairs within pairCutoff, sorted by distance.
        auto pairs = FilterResiduePairs(candidates, cfg.pairCutoff);

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
