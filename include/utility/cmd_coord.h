/**
 * cmd_coord.h
 *
 * The 'coord' subcommand: find all residues that have at least one atom within
 * a user-defined spherical region of space, specified by a centre coordinate
 * and a radius.
 *
 * This is useful when you already know the approximate location of a binding
 * site or region of interest in Cartesian space and want to enumerate every
 * residue present there across a set of PDB structures.
 *
 * Entry point: cactus::cmd::runCoord(argc, argv)
 *
 * Output format (grep-friendly):
 *   RESID  file=<path>  chain=<id>  resName=<name>  resSeq=<id>  minDist=<Å>
 */

#pragma once
#ifndef CMD_COORD_H
#define CMD_COORD_H

#include <string>
#include <vector>

namespace cactus::cmd {

// ── Config for the coord subcommand ──────────────────────────────────────────

struct CoordConfig {
    std::string inputPath  = "./";
    std::string outputFile = "coord_results.log";

    // Centre of the search sphere (Å). Must be set by the user.
    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;

    // Radius of the search sphere (Å).
    double radius = 10.0;

    // If set, only atoms belonging to these residue types are reported.
    // Empty means all residue types are accepted.
    std::vector<std::string> filterResidues = {};
};

/**
 * Run the coordinate-based search subcommand.
 *
 * @param argc  Argument count (argv[0] should be "coord").
 * @param argv  Argument vector.
 * @returns     0 on success, 1 on error.
 */
int runCoord(int argc, char* argv[]);

} // namespace cactus::cmd

#endif // CMD_COORD_H
