/**
 * main.cpp — subcommand dispatcher
 *
 * Usage:
 *   cactus_resid <subcommand> [OPTIONS]
 *
 * Subcommands:
 *   dist    Find pairs of acidic residues (ASP/GLU) near a carbohydrate
 *           by distance cutoff
 *   coord   Find residues within a spherical coordinate region
 *
 * Run `cactus_resid <subcommand> --help` for subcommand-specific options.
 */

#include "utility/cmd_dist.h"
#include "utility/cmd_coord.h"
#include <iostream>
#include <string>

static void printTopLevelHelp(const char* prog)
{
    std::cout
        << "Usage: " << prog << " <subcommand> [OPTIONS]\n\n"
        << "Subcommands:\n"
        << "  dist    Find pairs of acidic residues (ASP/GLU) near a carbohydrate\n"
        << "          by distance cutoff\n"
        << "  coord   Find residues within a spherical coordinate region\n\n"
        << "Run '" << prog << " <subcommand> --help' for subcommand-specific options.\n";
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        printTopLevelHelp(argv[0]);
        return 1;
    }

    const std::string subcmd = argv[1];

    // Shift argv so each subcommand sees its own name at argv[0]
    // and its flags starting at argv[1].
    if (subcmd == "dist")
        return cactus::cmd::runDist(argc - 1, argv + 1);

    if (subcmd == "coord")
        return cactus::cmd::runCoord(argc - 1, argv + 1);

    if (subcmd == "-h" || subcmd == "--help") {
        printTopLevelHelp(argv[0]);
        return 0;
    }

    std::cerr << "Unknown subcommand: '" << subcmd << "'\n"
              << "Run '" << argv[0] << " --help' for available subcommands.\n";
    return 1;
}
