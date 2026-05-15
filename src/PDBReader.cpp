#include <md/PDBReader.h>


using namespace mandacaru;

namespace mandacaru:md::pdb {

  struct Vec3 {
    double x, y, z;
  };

  struct PDBAtom {
      int         serial;      // atom serial number (1-based, as in PDB)
      std::string name;        // atom name  e.g. "CA"
      std::string resName;     // residue name
      char        chainID;
      int         resSeq;
      Vec3        pos;         // coordinates in Angstrom (PDB native)
      double      occupancy;
      double      mass;        // filled in by assignMasses()
  };

  
  // ─── PDB parser ───────────────────────────────────────────────────────────────
  // Follows the PDB fixed-column standard (ATOM/HETATM records).
  //
  //  Columns (1-based, inclusive):
  //   1-6   Record type
  //   7-11  Serial
  //  13-16  Atom name
  //  18-20  Residue name
  //  22     Chain ID
  //  23-26  Residue sequence number
  //  31-38  X (Å)
  //  39-46  Y (Å)
  //  47-54  Z (Å)
  //  55-60  Occupancy
  //  61-66  TempFactor
  //  77-78  Element symbol (often missing in GROMACS output — we derive from name)
  
  static std::string fixedCol(const std::string& line, int start, int end) {
      // start/end are 1-based, inclusive
      if (static_cast<int>(line.size()) < start) return "";
      int s = start - 1;
      int len = std::min(end - start + 1, static_cast<int>(line.size()) - s);
      std::string col = line.substr(s, len);
      // Trim whitespace
      auto a = col.find_first_not_of(' ');
      auto b = col.find_last_not_of(' ');
      return (a == std::string::npos) ? "" : col.substr(a, b - a + 1);
  }

  std::vector<PDBAtom> parsePDB(const std::string& path) {
      std::ifstream f(path);
      if (!f) throw std::runtime_error("parsePDB: cannot open " + path);

      std::vector<PDBAtom> atoms;
      std::string line;

      while (std::getline(f, line)) {
          if (line.size() < 27) continue;
          std::string rec = fixedCol(line, 1, 6);
          if (rec != "ATOM" && rec != "HETATM") continue;

          PDBAtom a;
          try {
              a.serial  = std::stoi(fixedCol(line, 7,  11));
              a.name    =           fixedCol(line, 13, 16);
              a.resName =           fixedCol(line, 18, 20);
              a.chainID = (line.size() >= 22) ? line[21] : ' ';
              a.resSeq  = std::stoi(fixedCol(line, 23, 26));
              a.pos.x   = std::stod(fixedCol(line, 31, 38));
              a.pos.y   = std::stod(fixedCol(line, 39, 46));
              a.pos.z   = std::stod(fixedCol(line, 47, 54));
              // Occupancy optional
              std::string occ = fixedCol(line, 55, 60);
              a.occupancy = occ.empty() ? 1.0 : std::stod(occ);
          } catch (...) {
              continue; // skip malformed records
          }

          a.mass = elementMass(a.name);
          atoms.push_back(a);
      }

      if (atoms.empty())
          throw std::runtime_error("parsePDB: no ATOM/HETATM records found in " + path);

      return atoms;
  }
}
