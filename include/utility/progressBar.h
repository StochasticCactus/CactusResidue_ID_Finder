/**
 * progressBar.h
 *
 * Lightweight header-only progress bar for terminal output.
 *
 * Usage:
 *
 *   size_t total = files.size();
 *   for (size_t i = 0; i < total; ++i) {
 *       printProgress(i, total, files[i].path());
 *       // ... do work ...
 *   }
 *   printProgress(total, total);   // finishes the bar with a newline
 */

#pragma once
#ifndef progressBar
#define progressBar

#include <iostream>
#include <string>
#include <algorithm>

namespace cactus::util {

/**
 * Prints (or updates) a progress bar on a single terminal line.
 *
 * Call once per iteration, then call a final time with `current == total`
 * to flush the completed bar and move to the next line.
 *
 * @param current   Number of items completed so far.
 * @param total     Total number of items.
 * @param label     Optional label shown after the bar (e.g. current filename).
 *                  Truncated to `labelWidth` characters so the line stays tidy.
 * @param barWidth  Visual width of the filled/empty bar in characters.
 * @param labelWidth Maximum characters reserved for the label column.
 */
inline void printProgress(size_t      current,
                          size_t      total,
                          std::string label      = "",
                          int         barWidth   = 40,
                          int         labelWidth = 30)
{
    if (total == 0) return;

    const double fraction = static_cast<double>(current) / static_cast<double>(total);
    const int    filled   = static_cast<int>(fraction * barWidth);
    const int    percent  = static_cast<int>(fraction * 100.0);

    // Truncate label so the line doesn't wrap on narrow terminals.
    if (static_cast<int>(label.size()) > labelWidth)
        label = "..." + label.substr(label.size() - (labelWidth - 3));

    std::cerr << '\r'
              << '['
              << std::string(filled,            '#')
              << std::string(barWidth - filled, '-')
              << ']'
              << ' ' << std::to_string(percent) << '%'
              << ' ' << std::to_string(current) << '/' << std::to_string(total)
              << "  " << std::left;

    // Pad label to labelWidth so leftover characters from a longer previous
    // label are overwritten when the filename gets shorter.
    std::string paddedLabel = label;
    paddedLabel.resize(labelWidth, ' ');
    std::cerr << paddedLabel << std::flush;

    // On the final call, end the line so subsequent output starts cleanly.
    if (current == total)
        std::cerr << '\n';
}

} // namespace cactus::util

#endif // progressBar
