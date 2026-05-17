/**
 * cmd_coord.h
 *
 * The 'coord' subcommand: find pairs of acidic residues (ASP / GLU) whose
 * carboxyl groups lie within a configurable distance of atoms in a reference
 * PDB structure, then report pairs of those candidates that are close to each
 * other.
 *
 * This is the same two-stage pipeline as 'dist', but instead of using
 * carbohydrate atoms as the reference set it uses a user-supplied reference
 * PDB file.  This is useful when you already have a structure that defines the
 * region of interest (e.g. a known inhibitor, a co-crystallised ligand, or
 * another protein chain) and want to find acidic residues that are
 * geometrically equivalent across a set of query structures.
 *
 * Entry point: cactus::cmd::runCoord(argc, argv)
 *
 * Usage:
 *   cactus_resid coord --ref reference.pdb [OPTIONS] [PATH]
 *
 * Output format (grep-friendly, same tag as 'dist' for unified downstream):
 *   PAIR  file=<path>  res1=<seqID>  res2=<seqID>  dist=<Å>
 */

#pragma once
#ifndef CMD_COORD_H
#define CMD_COORD_H

#include <string>
#include <vector>

namespace cactus::cmd {

// ── Config for the coord subcommand ──────────────────────────────────────────

/**
 * Run the reference-structure coordinate search subcommand.
 *
 * @param argc  Argument count (argv[0] should be "coord").
 * @param argv  Argument vector.
 * @returns     0 on success, 1 on error.
 */
int runCoord(int argc, char* argv[]);

} // namespace cactus::cmd

#endif // CMD_COORD_H
