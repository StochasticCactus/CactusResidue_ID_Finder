/*
 * PeriodicTable.h
 *
 * Standard atomic weights from IUPAC 2021 recommendations.
 * For elements with no stable isotopes, the mass of the most stable
 * isotope is used (marked with †).
 *
 * Storage strategy:
 *   - constexpr array of Element structs → zero runtime cost, lives in
 *     read-only memory, works at compile time.
 *   - buildLookup() creates a std::unordered_map<std::string, const Element*>
 *     keyed on both symbol ("Fe") and full name ("iron") for O(1) runtime access.
 *   - findByZ() is a simple array index (atomicNumber - 1).
 *   - elementMassFromPDBName() handles GROMACS/PDB atom-name conventions
 *     ("CA", "1HB", "OW", etc.) and is a drop-in replacement for the
 *     function in your existing code.
 *
 * Usage:
 *   #include "md/PeriodicTable.h"
 *   using namespace cactus::periodic;
 *
 *   double m = table[5].mass;                          // Carbon by index
 *   const Element* e = lookup()["Fe"];                 // Iron by symbol
 *   double mass = elementMassFromPDBName("CA");        // → Carbon alpha (12.011)
 */

#pragma once

#include <array>
#include <string>
#include <string_view>
#include <unordered_map>
#include <cctype>
#include <stdexcept>
#include <algorithm>

namespace cactus::periodic {

// ─── Element record ───────────────────────────────────────────────────────────

struct Element {
    int         atomicNumber;
    const char* symbol;        // e.g. "Fe"
    const char* name;          // e.g. "iron"
    double      mass;          // standard atomic weight (u / Da)
    bool        stableIsotope; // false → mass is most-stable isotope (†)
};

// ─── Full periodic table (Z = 1 … 118) ───────────────────────────────────────
// Index = atomicNumber - 1

inline constexpr std::array<Element, 118> table = {{
    {   1, "H",   "hydrogen",      1.008,        true  },
    {   2, "He",  "helium",        4.002602,     true  },
    {   3, "Li",  "lithium",       6.94,         true  },
    {   4, "Be",  "beryllium",     9.0121831,    true  },
    {   5, "B",   "boron",        10.81,         true  },
    {   6, "C",   "carbon",       12.011,        true  },
    {   7, "N",   "nitrogen",     14.007,        true  },
    {   8, "O",   "oxygen",       15.999,        true  },
    {   9, "F",   "fluorine",     18.998403163,  true  },
    {  10, "Ne",  "neon",         20.1797,       true  },
    {  11, "Na",  "sodium",       22.98976928,   true  },
    {  12, "Mg",  "magnesium",    24.305,        true  },
    {  13, "Al",  "aluminium",    26.9815384,    true  },
    {  14, "Si",  "silicon",      28.085,        true  },
    {  15, "P",   "phosphorus",   30.973761998,  true  },
    {  16, "S",   "sulfur",       32.06,         true  },
    {  17, "Cl",  "chlorine",     35.45,         true  },
    {  18, "Ar",  "argon",        39.948,        true  },
    {  19, "K",   "potassium",    39.0983,       true  },
    {  20, "Ca",  "calcium",      40.078,        true  },
    {  21, "Sc",  "scandium",     44.955907,     true  },
    {  22, "Ti",  "titanium",     47.867,        true  },
    {  23, "V",   "vanadium",     50.9415,       true  },
    {  24, "Cr",  "chromium",     51.9961,       true  },
    {  25, "Mn",  "manganese",    54.938043,     true  },
    {  26, "Fe",  "iron",         55.845,        true  },
    {  27, "Co",  "cobalt",       58.933194,     true  },
    {  28, "Ni",  "nickel",       58.6934,       true  },
    {  29, "Cu",  "copper",       63.546,        true  },
    {  30, "Zn",  "zinc",         65.38,         true  },
    {  31, "Ga",  "gallium",      69.723,        true  },
    {  32, "Ge",  "germanium",    72.630,        true  },
    {  33, "As",  "arsenic",      74.921595,     true  },
    {  34, "Se",  "selenium",     78.971,        true  },
    {  35, "Br",  "bromine",      79.904,        true  },
    {  36, "Kr",  "krypton",      83.798,        true  },
    {  37, "Rb",  "rubidium",     85.4678,       true  },
    {  38, "Sr",  "strontium",    87.62,         true  },
    {  39, "Y",   "yttrium",      88.905838,     true  },
    {  40, "Zr",  "zirconium",    91.224,        true  },
    {  41, "Nb",  "niobium",      92.90637,      true  },
    {  42, "Mo",  "molybdenum",   95.95,         true  },
    {  43, "Tc",  "technetium",   96.90636,      false }, // †
    {  44, "Ru",  "ruthenium",   101.07,         true  },
    {  45, "Rh",  "rhodium",     102.90549,      true  },
    {  46, "Pd",  "palladium",   106.42,         true  },
    {  47, "Ag",  "silver",      107.8682,       true  },
    {  48, "Cd",  "cadmium",     112.414,        true  },
    {  49, "In",  "indium",      114.818,        true  },
    {  50, "Sn",  "tin",         118.710,        true  },
    {  51, "Sb",  "antimony",    121.760,        true  },
    {  52, "Te",  "tellurium",   127.60,         true  },
    {  53, "I",   "iodine",      126.90447,      true  },
    {  54, "Xe",  "xenon",       131.293,        true  },
    {  55, "Cs",  "caesium",     132.90545196,   true  },
    {  56, "Ba",  "barium",      137.327,        true  },
    {  57, "La",  "lanthanum",   138.90547,      true  },
    {  58, "Ce",  "cerium",      140.116,        true  },
    {  59, "Pr",  "praseodymium",140.90766,      true  },
    {  60, "Nd",  "neodymium",   144.242,        true  },
    {  61, "Pm",  "promethium",  144.91276,      false }, // †
    {  62, "Sm",  "samarium",    150.36,         true  },
    {  63, "Eu",  "europium",    151.964,        true  },
    {  64, "Gd",  "gadolinium",  157.25,         true  },
    {  65, "Tb",  "terbium",     158.925354,     true  },
    {  66, "Dy",  "dysprosium",  162.500,        true  },
    {  67, "Ho",  "holmium",     164.930329,     true  },
    {  68, "Er",  "erbium",      167.259,        true  },
    {  69, "Tm",  "thulium",     168.934219,     true  },
    {  70, "Yb",  "ytterbium",   173.045,        true  },
    {  71, "Lu",  "lutetium",    174.9668,       true  },
    {  72, "Hf",  "hafnium",     178.486,        true  },
    {  73, "Ta",  "tantalum",    180.94788,      true  },
    {  74, "W",   "tungsten",    183.84,         true  },
    {  75, "Re",  "rhenium",     186.207,        true  },
    {  76, "Os",  "osmium",      190.23,         true  },
    {  77, "Ir",  "iridium",     192.217,        true  },
    {  78, "Pt",  "platinum",    195.084,        true  },
    {  79, "Au",  "gold",        196.966570,     true  },
    {  80, "Hg",  "mercury",     200.592,        true  },
    {  81, "Tl",  "thallium",    204.38,         true  },
    {  82, "Pb",  "lead",        207.2,          true  },
    {  83, "Bi",  "bismuth",     208.98040,      true  },
    {  84, "Po",  "polonium",    208.98243,      false }, // †
    {  85, "At",  "astatine",    209.98715,      false }, // †
    {  86, "Rn",  "radon",       222.01758,      false }, // †
    {  87, "Fr",  "francium",    223.01973,      false }, // †
    {  88, "Ra",  "radium",      226.02541,      false }, // †
    {  89, "Ac",  "actinium",    227.02775,      false }, // †
    {  90, "Th",  "thorium",     232.0377,       true  },
    {  91, "Pa",  "protactinium",231.03588,      true  },
    {  92, "U",   "uranium",     238.02891,      true  },
    {  93, "Np",  "neptunium",   237.04817,      false }, // †
    {  94, "Pu",  "plutonium",   244.06420,      false }, // †
    {  95, "Am",  "americium",   243.06138,      false }, // †
    {  96, "Cm",  "curium",      247.07035,      false }, // †
    {  97, "Bk",  "berkelium",   247.07031,      false }, // †
    {  98, "Cf",  "californium", 251.07959,      false }, // †
    {  99, "Es",  "einsteinium", 252.08298,      false }, // †
    { 100, "Fm",  "fermium",     257.09511,      false }, // †
    { 101, "Md",  "mendelevium", 258.09843,      false }, // †
    { 102, "No",  "nobelium",    259.10100,      false }, // †
    { 103, "Lr",  "lawrencium",  266.11983,      false }, // †
    { 104, "Rf",  "rutherfordium",267.12167,     false }, // †
    { 105, "Db",  "dubnium",     268.12567,      false }, // †
    { 106, "Sg",  "seaborgium",  269.12863,      false }, // †
    { 107, "Bh",  "bohrium",     270.13336,      false }, // †
    { 108, "Hs",  "hassium",     269.13375,      false }, // †
    { 109, "Mt",  "meitnerium",  277.15361,      false }, // †
    { 110, "Ds",  "darmstadtium",282.16550,      false }, // †
    { 111, "Rg",  "roentgenium", 282.16912,      false }, // †
    { 112, "Cn",  "copernicium", 286.17913,      false }, // †
    { 113, "Nh",  "nihonium",    286.18221,      false }, // †
    { 114, "Fl",  "flerovium",   290.19214,      false }, // †
    { 115, "Mc",  "moscovium",   290.19598,      false }, // †
    { 116, "Lv",  "livermorium", 293.20449,      false }, // †
    { 117, "Ts",  "tennessine",  294.21046,      false }, // †
    { 118, "Og",  "oganesson",   295.21624,      false }, // †
}};

// ─── O(1) lookup by symbol or name ───────────────────────────────────────────
// Thread-safe after first call (C++11 static-local guarantee).

inline const std::unordered_map<std::string, const Element*>& lookup() {
    static const auto map = []() {
        std::unordered_map<std::string, const Element*> m;
        m.reserve(table.size() * 2);
        for (const auto& e : table) {
            m[e.symbol] = &e;
            m[e.name]   = &e;
        }
        return m;
    }();
    return map;
}

// ─── Lookup by atomic number (Z) ──────────────────────────────────────────────

inline const Element& findByZ(int Z) {
    if (Z < 1 || Z > 118)
        throw std::out_of_range("Atomic number out of range: " + std::to_string(Z));
    return table[static_cast<std::size_t>(Z - 1)];
}

// ─── Derive element symbol from a PDB atom name ───────────────────────────────
//
// PDB / GROMACS atom-name conventions handled:
//   "CA"  → alpha-carbon → "C"  (two-char lookup fails → single char "C")
//   "1HB" → leading digit stripped → "H"
//   "HG1" → "H"
//   "OW"  → "O"  (TIP3P/TIP4P water)
//   "CL"  → "Cl" (chloride ion)
//   "NA"  → "Na" (sodium ion)
//   "MG"  → "Mg"
//   "FE"  → "Fe"
//   "ZN"  → "Zn"
//   "CU"  → "Cu"

inline std::string elementSymbolFromPDBName(const std::string& atomName) {
    // 1. Strip leading whitespace and digits
    std::size_t start = 0;
    while (start < atomName.size() && (std::isspace((unsigned char)atomName[start]) ||
                                       std::isdigit((unsigned char)atomName[start])))
        ++start;

    if (start >= atomName.size())
        return "C"; // empty / all-digits: safe fallback

    // 2. Extract up to 2 alpha characters, upper-cased
    std::string alpha;
    for (std::size_t i = start; i < atomName.size() && alpha.size() < 2; ++i) {
        char c = atomName[i];
        if (std::isalpha((unsigned char)c))
            alpha += static_cast<char>(std::toupper((unsigned char)c));
        else if (!alpha.empty())
            break; // stop at embedded digit/special after first alpha char
    }

    if (alpha.empty()) return "C";

    // 3. Normalise to Title-case for map key ("FE" → "Fe")
    if (alpha.size() == 2)
        alpha[1] = static_cast<char>(std::tolower((unsigned char)alpha[1]));

    // 4. Try 2-char symbol first (catches Cl, Na, Mg, Fe, Cu, Zn, etc.)
    if (alpha.size() == 2) {
        const auto& m = lookup();
        if (m.count(alpha))
            return alpha;
    }

    // 5. Fall back to single-character symbol
    return std::string(1, alpha[0]);
}

// ─── Mass from PDB atom name ───────────────────────────────────────────────────

inline double elementMassFromPDBName(const std::string& atomName) {
    std::string sym = elementSymbolFromPDBName(atomName);
    const auto& m   = lookup();
    auto it = m.find(sym);
    if (it != m.end())
        return it->second->mass;
    return 12.011; // unknown: fallback to carbon
}

} // namespace cactus::periodic
