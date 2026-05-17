/**
 * cmd_dist.h
 *
 * The 'dist' subcommand: find pairs of acidic residues (ASP / GLU) whose
 * catalytic oxygens lie within a configurable distance of a carbohydrate
 * residue, then report pairs of those candidates that are close to each other.
 *
 * Entry point: cactus::cmd::runDist(argc, argv)
 */

#pragma once
#ifndef CMD_DIST_H
#define CMD_DIST_H

namespace cactus::cmd {

/**
 * Run the distance-based search subcommand.
 *
 * @param argc  Argument count (argv[0] should be "dist").
 * @param argv  Argument vector.
 * @returns     0 on success, 1 on error (suitable for returning from main).
 */
int runDist(int argc, char* argv[]);

} // namespace cactus::cmd

#endif // CMD_DIST_H
