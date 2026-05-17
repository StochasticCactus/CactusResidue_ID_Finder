/**
 * cmd_coord.cpp
 *
 * Implementation of the 'coord' subcommand.
 *
 * For each PDB file in the input directory, every residue that has at least
 * one atom within `radius` Å of the user-supplied centre (cx, cy, cz) is
 * reported.  Each hit is written as a single grep-able line:
 *
 *   RESID  file=<path>  chain=<id>  resName=<name>  resSeq=<id>  minDist=<Å>
 *
 * where minDist is the distance from the centre to the closest atom of that
 * residue.  Results within a file are sorted by ascending minDist.
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
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace cactus::pdb;
using namespace cactus::sort;

namespace cactus::cmd {

// ── Argument parsing ──────────────────────────────────────────────────────────

static void printCoordHelp(const char* prog)
{
    std::cout
        << "Usage: " << prog << " [OPTIONS] [PATH]\n\n"
        << "Find all residues with at least one atom within a spherical region\n"
        << "of space, across all PDB files in PATH.\n\n"
        << "Arguments:\n"
        << "  PATH                    Directory containing PDB files (default: ./)\n\n"
        << "Options:\n"
        << "  -h, --help              Show this help message and exit\n"
        << "  -o, --output FILE       Write results to FILE (default: coord_results.log)\n"
        << "  --cx X                  X coordinate of the sphere centre (default: 0.0)\n"
        << "  --cy Y                  Y coordinate of the sphere centre (default: 0.0)\n"
        << "  --cz Z                  Z coordinate of the sphere centre (default: 0.0)\n"
        << "  --radius R              Search sphere radius in Å (default: 10.0)\n"
        << "  --residues RES          Comma-separated residue names to restrict output to\n"
        << "                          (default: all residues)\n\n"
        << "Output format:\n"
        << "  RESID  file=<path>  chain=<id>  resName=<name>  resSeq=<id>  minDist=<Å>\n\n"
        << "  Recover all hits with:\n"
        << "    grep '^RESID' coord_results.log\n\n"
        << "Examples:\n"
        << "  " << prog << " ./structures/ --cx 12.3 --cy 45.6 --cz 7.8 --radius 8.0\n"
        << "  " << prog << " ./structures/ --cx 0 --cy 0 --cz 0 --radius 15 --residues GLU,ASP\n";
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

        } else if (arg == "--cx" && i + 1 < argc) {
            cfg.cx = std::stod(argv[++i]);

        } else if (arg == "--cy" && i + 1 < argc) {
            cfg.cy = std::stod(argv[++i]);

        } else if (arg == "--cz" && i + 1 < argc) {
            cfg.cz = std::stod(argv[++i]);

        } else if (arg == "--radius" && i + 1 < argc) {
            cfg.radius = std::stod(argv[++i]);

        } else if (arg == "--residues" && i + 1 < argc) {
            cfg.filterResidues = splitComma(argv[++i]);

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

// ── Search logic ──────────────────────────────────────────────────────────────

// A residue hit: the closest atom distance and enough identity fields to
// write a useful output line.
struct ResidHit {
    char        chainID;
    std::string resName;
    int         resSeq;
    double      minDist;  // distance from sphere centre to closest atom
};

/**
 * For a given atom list, find every unique residue that has at least one atom
 * within `radius` of the point (cx, cy, cz).
 *
 * Keyed on {chainID, resSeq} so residues with the same sequence number on
 * different chains are treated as distinct.  Only the minimum per-atom
 * distance is stored for each residue.
 */
static std::vector<ResidHit> findResiduesInSphere(
        const std::vector<PDBAtom>& atoms,
        double cx, double cy, double cz, double radius,
        const std::vector<std::string>& filterResidues)
{
    const Vec3 centre{cx, cy, cz};

    // Key: {chainID, resSeq}  →  best (minimum) distance seen for that residue.
    std::map<std::pair<char, int>, ResidHit> best;

    for (const auto& atom : atoms) {

        // If the user asked for specific residue types, skip others.
        if (!filterResidues.empty()) {
            bool wanted = false;
            for (const auto& r : filterResidues)
                if (atom.resName == r) { wanted = true; break; }
            if (!wanted) continue;
        }

        const double dist = distance(centre, atom.pos);
        if (dist >= radius) continue;

        const auto key = std::make_pair(atom.chainID, atom.resSeq);
        auto it = best.find(key);
        if (it == best.end() || dist < it->second.minDist)
            best[key] = {atom.chainID, atom.resName, atom.resSeq, dist};
    }

    // Flatten and sort by ascending distance.
    std::vector<ResidHit> hits;
    hits.reserve(best.size());
    for (const auto& kv : best)
        hits.push_back(kv.second);

    std::sort(hits.begin(), hits.end(),
        [](const ResidHit& a, const ResidHit& b) {
            return a.minDist < b.minDist;
        });

    return hits;
}

// ── Entry point ───────────────────────────────────────────────────────────────

int runCoord(int argc, char* argv[])
{
    const CoordConfig cfg = parseCoordArgs(argc, argv);

    std::ofstream logFile(cfg.outputFile);
    if (!logFile) {
        std::cerr << "Error: could not open output file: " << cfg.outputFile << "\n";
        return 1;
    }

    // Self-describing log header.
    logFile << "# [coord] Residues within a spherical coordinate region\n"
            << "# centre : (" << cfg.cx << ", " << cfg.cy << ", " << cfg.cz << ") Å\n"
            << "# radius : " << cfg.radius << " Å\n";
    if (!cfg.filterResidues.empty()) {
        logFile << "# filter : ";
        for (size_t k = 0; k < cfg.filterResidues.size(); ++k) {
            logFile << cfg.filterResidues[k];
            if (k + 1 < cfg.filterResidues.size()) logFile << ", ";
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

        const auto hits = findResiduesInSphere(
            atoms, cfg.cx, cfg.cy, cfg.cz, cfg.radius, cfg.filterResidues);

        logFile << "# FILE: " << file.path() << "\n";

        if (hits.empty()) {
            logFile << "  (no residues found)\n";
            continue;
        }

        for (const auto& h : hits) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(3);
            oss << "RESID"
                << "  file="    << file.path()
                << "  chain="   << h.chainID
                << "  resName=" << h.resName
                << "  resSeq="  << h.resSeq
                << "  minDist=" << h.minDist;
            logFile << oss.str() << "\n";
        }
    }

    cactus::util::printProgress(total, total);
    std::cout << "Results written to: " << cfg.outputFile << "\n";
    return 0;
}

} // namespace cactus::cmd
