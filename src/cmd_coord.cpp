/**
 * cmd_coord.cpp
 *
 * Implementation of the 'coord' subcommand.
 *
 * Given a set of target points in space, this tool searches every query PDB
 * for the single acidic residue (ASP or GLU) whose carboxyl-group midpoint is
 * nearest to each target point.  One MATCH or NOMATCH line is written per
 * target point per file, making it straightforward to recover residue numbers
 * across large model sets (e.g. RFDiffusion outputs) where catalytic residues
 * are geometrically fixed but their sequence numbers vary.
 *
 * Target points can be supplied in two non-exclusive ways:
 *
 *   --ref FILE --ref-seqid 142,158
 *       Parse a reference PDB and compute the carboxyl midpoint of each
 *       listed residue.  Works for both ASP (OD1/OD2) and GLU (OE1/OE2).
 *
 *   --point x,y,z   (repeatable)
 *       Supply an explicit coordinate directly on the command line.
 *
 * Both can be combined; all points are pooled and numbered in the order they
 * are collected (reference residues first, then --point arguments).
 *
 * Matching logic
 * ──────────────
 * For each target point the tool finds the single ASP/GLU carboxyl midpoint
 * that is closest.  If that minimum distance exceeds --cutoff the result is
 * flagged NOMATCH; otherwise MATCH is written with the residue's sequence
 * number, residue name, and distance.
 *
 * Output (grep-friendly)
 * ──────────────────────
 *   MATCH    file=<path>  point=<N>  res=<seqID>  resname=<ASP|GLU>  dist=<Å>
 *   NOMATCH  file=<path>  point=<N>  nearest=<Å>  cutoff=<Å>
 */

#include "utility/cmd_coord.h"
#include "utility/fileHandler.h"
#include "utility/progressBar.h"
#include "md/PDBReader.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cactus::pdb;

namespace cactus::cmd {

// ─────────────────────────────────────────────────────────────────────────────
// Internal types
// ─────────────────────────────────────────────────────────────────────────────

/** A target point together with a human-readable label for log output. */
struct TargetPoint {
    Vec3        coords;
    std::string label;   // e.g. "ref:ASP142" or "point:5.0,2.0,7.0"
};

/** Everything that describes one candidate acidic residue in a query file. */
struct AcidicMidpoint {
    int         resSeq;
    std::string resName;
    Vec3        midpoint;
};

// ─────────────────────────────────────────────────────────────────────────────
// Geometry helpers
// ─────────────────────────────────────────────────────────────────────────────

static double dist3(const Vec3& a, const Vec3& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// ─────────────────────────────────────────────────────────────────────────────
// Carboxyl oxygen names for each supported acidic residue
// ─────────────────────────────────────────────────────────────────────────────

static const std::unordered_map<std::string, std::array<std::string, 2>>
CARBOXYL_OXYGENS = {
    {"ASP", {"OD1", "OD2"}},
    {"GLU", {"OE1", "OE2"}}
};

// ─────────────────────────────────────────────────────────────────────────────
// Argument parsing
// ─────────────────────────────────────────────────────────────────────────────

static void printCoordHelp(const char* prog)
{
    std::cout
        << "Usage: " << prog << " coord [OPTIONS] [PATH]\n\n"
        << "For each query PDB in PATH, find the single ASP/GLU residue whose\n"
        << "carboxyl-group midpoint is nearest to each supplied target point.\n"
        << "Designed for locating catalytic residues across large model sets\n"
        << "(e.g. RFDiffusion outputs) where coordinates are fixed but residue\n"
        << "numbers vary.\n\n"
        << "Target points (at least one required; both sources may be combined):\n"
        << "  --ref FILE              Reference PDB from which to derive target\n"
        << "                          points (use with --ref-seqid)\n"
        << "  --ref-seqid IDS         Comma-separated residue sequence numbers in\n"
        << "                          the reference whose carboxyl midpoints become\n"
        << "                          target points (e.g. 142,158)\n"
        << "  --point x,y,z           Explicit target coordinate in Å (repeatable)\n\n"
        << "Options:\n"
        << "  -h, --help              Show this help message and exit\n"
        << "  -o, --output FILE       Write results to FILE\n"
        << "                          (default: coord_results.log)\n"
        << "  --cutoff DIST           Maximum distance in Å for a valid match;\n"
        << "                          farther results are flagged NOMATCH\n"
        << "                          (default: 5.0)\n"
        << "  --residues RES          Comma-separated acidic residue types to\n"
        << "                          search in query structures (default: GLU,ASP)\n\n"
        << "Arguments:\n"
        << "  PATH                    Directory of query PDB files (default: ./)\n\n"
        << "Output format (one line per target point per file):\n"
        << "  MATCH    file=<path>  point=<N>  res=<seqID>  resname=<name>  dist=<Å>\n"
        << "  NOMATCH  file=<path>  point=<N>  nearest=<Å>  cutoff=<Å>\n\n"
        << "Grep recipes:\n"
        << "  All successful matches:          grep '^MATCH'  coord_results.log\n"
        << "  All failures:                    grep '^NOMATCH' coord_results.log\n"
        << "  Matches for point 2 only:        grep '^MATCH.*point=2' coord_results.log\n"
        << "  Files where any point failed:    grep '^NOMATCH' coord_results.log | \\\n"
        << "                                     awk -F'file=' '{print $2}' | cut -d' ' -f1\n\n"
        << "Examples:\n"
        << "  # Derive two target points from a reference PDB:\n"
        << "  " << prog << " coord --ref ref.pdb --ref-seqid 142,158 ./models/\n\n"
        << "  # Supply target coordinates explicitly:\n"
        << "  " << prog << " coord --point 5.0,2.0,7.0 --point 5.0,3.0,8.0 ./models/\n\n"
        << "  # Mix both sources, wider cutoff, GLU only:\n"
        << "  " << prog << " coord --ref ref.pdb --ref-seqid 142 \\\n"
        << "          --point 5.0,3.0,8.0 --cutoff 8.0 --residues GLU ./models/\n";
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
    for (const auto& tok : splitComma(s))
        result.push_back(std::stoi(tok));
    return result;
}

/** Parse "--point x,y,z" into a Vec3.  Throws std::invalid_argument on bad input. */
static Vec3 parsePoint(const std::string& s)
{
    const auto parts = splitComma(s);
    if (parts.size() != 3)
        throw std::invalid_argument("--point requires exactly three comma-separated values: " + s);
    return { std::stod(parts[0]), std::stod(parts[1]), std::stod(parts[2]) };
}

struct CoordConfig {
    std::string              inputPath   = "./";
    std::string              outputFile  = "coord_results.log";
    std::string              refPath;
    std::vector<int>         refSeqIDs;
    std::vector<Vec3>        explicitPoints;   // from --point
    double                   cutoff      = 5.0;
    std::vector<std::string> acidicResidues = {"GLU", "ASP"};
};

static CoordConfig parseCoordArgs(int argc, char* argv[])
{
    CoordConfig cfg;

    for (int i = 2; i < argc; ++i) {   // argv[1] is the subcommand "coord"
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printCoordHelp(argv[0]);
            std::exit(0);

        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            cfg.outputFile = argv[++i];

        } else if (arg == "--ref" && i + 1 < argc) {
            cfg.refPath = argv[++i];

        } else if (arg == "--ref-seqid" && i + 1 < argc) {
            cfg.refSeqIDs = splitInt(argv[++i]);

        } else if (arg == "--point" && i + 1 < argc) {
            cfg.explicitPoints.push_back(parsePoint(argv[++i]));

        } else if (arg == "--cutoff" && i + 1 < argc) {
            cfg.cutoff = std::stod(argv[++i]);

        } else if (arg == "--residues" && i + 1 < argc) {
            cfg.acidicResidues = splitComma(argv[++i]);

        } else if (arg[0] != '-') {
            cfg.inputPath = arg;

        } else {
            std::cerr << "Unknown option: " << arg
                      << "  (run '" << argv[0] << " coord --help' for usage)\n";
            std::exit(1);
        }
    }
    return cfg;
}

// ─────────────────────────────────────────────────────────────────────────────
// Derive target points from reference PDB
// ─────────────────────────────────────────────────────────────────────────────

/**
 * For each requested sequence ID, find the matching residue in the reference
 * PDB, compute its carboxyl midpoint, and return it as a TargetPoint.
 *
 * Only ASP (OD1/OD2) and GLU (OE1/OE2) are supported as reference residues.
 * If a requested seqID is not found, or the residue is not an acidic type with
 * known carboxyl oxygens, an error is printed and the program exits.
 */
static std::vector<TargetPoint> targetPointsFromRef(
        const std::string&      refPath,
        const std::vector<int>& seqIDs)
{
    const std::vector<PDBAtom> allAtoms = parsePDB(refPath);

    std::vector<TargetPoint> targets;

    for (const int id : seqIDs) {

        // Find the carboxyl oxygen names expected for this residue.
        // First pass: identify the residue name for this seqID.
        std::string resName;
        for (const auto& a : allAtoms) {
            if (a.resSeq == id) { resName = a.resName; break; }
        }

        if (resName.empty()) {
            std::cerr << "Error: residue with seqID " << id
                      << " not found in reference PDB '" << refPath << "'.\n";
            std::exit(1);
        }

        const auto it = CARBOXYL_OXYGENS.find(resName);
        if (it == CARBOXYL_OXYGENS.end()) {
            std::cerr << "Error: residue " << resName << " (seqID " << id
                      << ") in the reference is not a supported acidic residue "
                      << "(ASP or GLU).\n";
            std::exit(1);
        }

        const std::string& o1name = it->second[0];
        const std::string& o2name = it->second[1];

        // Collect the two carboxyl oxygens.
        Vec3 o1{}, o2{};
        bool found1 = false, found2 = false;

        for (const auto& a : allAtoms) {
            if (a.resSeq != id) continue;
            if (a.name == o1name) { o1 = {a.pos.x, a.pos.y, a.pos.z}; found1 = true; }
            if (a.name == o2name) { o2 = {a.pos.x, a.pos.y, a.pos.z}; found2 = true; }
            if (found1 && found2) break;
        }

        if (!found1 || !found2) {
            std::cerr << "Error: could not find both carboxyl oxygens ("
                      << o1name << ", " << o2name << ") for " << resName
                      << " seqID " << id << " in '" << refPath << "'.\n";
            std::exit(1);
        }

        const Vec3 midpoint = {
            (o1.x + o2.x) * 0.5,
            (o1.y + o2.y) * 0.5,
            (o1.z + o2.z) * 0.5
        };

        std::ostringstream label;
        label << "ref:" << resName << id;

        targets.push_back({ midpoint, label.str() });
    }

    return targets;
}

// ─────────────────────────────────────────────────────────────────────────────
// Extract acidic midpoints from a query PDB
// ─────────────────────────────────────────────────────────────────────────────

/**
 * For every residue in 'atoms' whose name is in 'allowedResidues', compute the
 * carboxyl midpoint and return it.  Residues that are missing one or both
 * carboxyl oxygens are silently skipped (malformed HETATM/ATOM records).
 */
static std::vector<AcidicMidpoint> extractAcidicMidpoints(
        const std::vector<PDBAtom>&      atoms,
        const std::vector<std::string>&  allowedResidues)
{
    // Group atoms by (resSeq, resName) so we can find both oxygens per residue.
    // Key: resSeq (sufficient — sequence numbers are unique within a chain for
    // the model files we target).
    struct ResData {
        std::string resName;
        std::unordered_map<std::string, Vec3> atomCoords;
    };
    std::unordered_map<int, ResData> bySeq;

    for (const auto& a : atoms) {
        // Filter to allowed residue types early.
        bool allowed = false;
        for (const auto& r : allowedResidues)
            if (a.resName == r) { allowed = true; break; }
        if (!allowed) continue;

        auto& rd = bySeq[a.resSeq];
        rd.resName = a.resName;
        rd.atomCoords[a.name] = { a.pos.x, a.pos.y, a.pos.z };
    }

    std::vector<AcidicMidpoint> result;

    for (const auto& [seq, rd] : bySeq) {
        const auto it = CARBOXYL_OXYGENS.find(rd.resName);
        if (it == CARBOXYL_OXYGENS.end()) continue;

        const std::string& o1name = it->second[0];
        const std::string& o2name = it->second[1];

        const auto a1 = rd.atomCoords.find(o1name);
        const auto a2 = rd.atomCoords.find(o2name);
        if (a1 == rd.atomCoords.end() || a2 == rd.atomCoords.end()) continue;

        const Vec3 mid = {
            (a1->second.x + a2->second.x) * 0.5,
            (a1->second.y + a2->second.y) * 0.5,
            (a1->second.z + a2->second.z) * 0.5
        };

        result.push_back({ seq, rd.resName, mid });
    }

    // Sort by sequence number for deterministic output.
    std::sort(result.begin(), result.end(),
              [](const AcidicMidpoint& a, const AcidicMidpoint& b) {
                  return a.resSeq < b.resSeq;
              });

    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Entry point
// ─────────────────────────────────────────────────────────────────────────────

int runCoord(int argc, char* argv[])
{
    const CoordConfig cfg = parseCoordArgs(argc, argv);

    // ── Validate acidic residue list ─────────────────────────────────────────

    for (const auto& r : cfg.acidicResidues) {
        if (CARBOXYL_OXYGENS.find(r) == CARBOXYL_OXYGENS.end()) {
            std::cerr << "Warning: '" << r
                      << "' is not a known acidic residue (ASP/GLU) — skipping.\n";
        }
    }

    // ── Build target point list ───────────────────────────────────────────────

    std::vector<TargetPoint> targets;

    // 1. From reference PDB + seqIDs.
    if (!cfg.refPath.empty()) {
        if (cfg.refSeqIDs.empty()) {
            std::cerr << "Error: --ref was given but --ref-seqid is missing.\n"
                      << "Specify which residues to use, e.g. --ref-seqid 142,158\n";
            return 1;
        }
        try {
            const auto refTargets = targetPointsFromRef(cfg.refPath, cfg.refSeqIDs);
            targets.insert(targets.end(), refTargets.begin(), refTargets.end());
        } catch (const std::exception& ex) {
            std::cerr << "Error reading reference file '" << cfg.refPath
                      << "': " << ex.what() << "\n";
            return 1;
        }
    }

    // 2. From explicit --point arguments.
    for (size_t i = 0; i < cfg.explicitPoints.size(); ++i) {
        const Vec3& p = cfg.explicitPoints[i];
        std::ostringstream label;
        label << std::fixed << std::setprecision(3)
              << "point:" << p.x << "," << p.y << "," << p.z;
        targets.push_back({ p, label.str() });
    }

    if (targets.empty()) {
        std::cerr << "Error: no target points defined.\n"
                  << "Use --ref + --ref-seqid and/or --point x,y,z.\n"
                  << "Run '" << argv[0] << " coord --help' for usage.\n";
        return 1;
    }

    std::cout << targets.size() << " target point(s):\n";
    for (size_t i = 0; i < targets.size(); ++i) {
        const auto& t = targets[i];
        std::cout << "  point " << (i + 1) << "  (" << t.label << ")  "
                  << std::fixed << std::setprecision(3)
                  << t.coords.x << "  " << t.coords.y << "  " << t.coords.z << "\n";
    }

    // ── Open log file ─────────────────────────────────────────────────────────

    std::ofstream logFile(cfg.outputFile);
    if (!logFile) {
        std::cerr << "Error: could not open output file: " << cfg.outputFile << "\n";
        return 1;
    }

    // Self-describing header.
    logFile << "# [coord] Nearest acidic residue to each target point\n"
            << "# cutoff    : " << cfg.cutoff << " Å\n"
            << "# residues  : ";
    for (size_t k = 0; k < cfg.acidicResidues.size(); ++k) {
        logFile << cfg.acidicResidues[k];
        if (k + 1 < cfg.acidicResidues.size()) logFile << ", ";
    }
    logFile << "\n# points    :\n";
    for (size_t i = 0; i < targets.size(); ++i) {
        const auto& t = targets[i];
        logFile << "#   " << (i + 1) << "  " << t.label << "  "
                << std::fixed << std::setprecision(3)
                << t.coords.x << "  " << t.coords.y << "  " << t.coords.z << "\n";
    }
    logFile << "#\n";

    // ── Query loop ────────────────────────────────────────────────────────────

    StructuresDir wdir(cfg.inputPath);
    const std::vector<std::filesystem::directory_entry> files(
        wdir.dirFiles(), std::filesystem::directory_iterator{});
    const size_t total = files.size();

    size_t totalMatches  = 0;
    size_t totalNoMatch  = 0;

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

        const auto midpoints = extractAcidicMidpoints(atoms, cfg.acidicResidues);

        logFile << "# FILE: " << file.path() << "\n";

        if (midpoints.empty()) {
            // No acidic residues at all — every point is a NOMATCH.
            for (size_t pi = 0; pi < targets.size(); ++pi) {
                logFile << "NOMATCH"
                        << "  file=" << file.path()
                        << "  point=" << (pi + 1)
                        << "  nearest=N/A"
                        << "  cutoff=" << std::fixed << std::setprecision(3)
                        << cfg.cutoff << "\n";
                ++totalNoMatch;
            }
            continue;
        }

        // For each target point: find the nearest acidic midpoint.
        for (size_t pi = 0; pi < targets.size(); ++pi) {
            const Vec3& tgt = targets[pi].coords;

            double      bestDist = std::numeric_limits<double>::max();
            const AcidicMidpoint* bestRes = nullptr;

            for (const auto& mp : midpoints) {
                const double d = dist3(tgt, mp.midpoint);
                if (d < bestDist) {
                    bestDist = d;
                    bestRes  = &mp;
                }
            }

            std::ostringstream line;
            line << std::fixed << std::setprecision(3);

            if (bestDist <= cfg.cutoff) {
                line << "MATCH"
                     << "  file="    << file.path()
                     << "  point="   << (pi + 1)
                     << "  res="     << bestRes->resSeq
                     << "  resname=" << bestRes->resName
                     << "  dist="    << bestDist;
                ++totalMatches;
            } else {
                line << "NOMATCH"
                     << "  file="    << file.path()
                     << "  point="   << (pi + 1)
                     << "  nearest=" << bestDist
                     << "  cutoff="  << cfg.cutoff;
                ++totalNoMatch;
            }

            logFile << line.str() << "\n";
        }
    }

    cactus::util::printProgress(total, total);

    std::cout << "\nDone.  " << totalMatches << " match(es), "
              << totalNoMatch << " NOMATCH(es) across "
              << total << " file(s).\n"
              << "Results written to: " << cfg.outputFile << "\n";

    return 0;
}

} // namespace cactus::cmd
