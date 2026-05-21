#pragma once
#ifndef CMD_COORD_H
#define CMD_COORD_H

namespace cactus::cmd {

/**
 * 'coord' subcommand entry point.
 *
 * For each query PDB in a directory, finds the single ASP/GLU residue whose
 * carboxyl-group midpoint is nearest to each supplied target point.  Target
 * points are derived from a reference PDB (--ref + --ref-seqid), from explicit
 * coordinates (--point x,y,z), or both combined.
 *
 * One MATCH or NOMATCH line is written per target point per file.
 *
 * @param argc  Argument count (argv[0] == "cactus_resid", argv[1] == "coord").
 * @param argv  Argument vector.
 * @returns     0 on success, 1 on error.
 */
int runCoord(int argc, char* argv[]);

} // namespace cactus::cmd

#endif // CMD_COORD_H
