# cactus_resid

A command-line tool for finding pairs of acidic amino acid residues (aspartate and glutamate) that are geometrically positioned near a region of interest in protein structures. It works on collections of PDB files — the standard file format used to store 3D protein structures.

---

## What problem does this solve?

In many enzymes, pairs of acidic residues (aspartate / ASP and glutamate / GLU) work together as a catalytic team. Their side chains each end in a carboxyl group (two oxygens), and the geometry of that pair — how close they are, and where they sit in the structure — is often critical to the enzyme's function.

This tool automates the process of finding such pairs across many structures at once, and can check whether those pairs are located near a specific chemical feature (a sugar, a ligand, another residue, etc.).

---

## Installation

The tool is written in C++ and needs to be compiled before use. If your collaborator has already sent you a compiled binary, skip this section entirely and go straight to [Usage](#usage).

If you need to compile it yourself you will need:

- A C++ compiler (GCC 11+ or Clang 13+)
- CMake 3.16 or newer

```bash
git clone <repository-url>
cd CactusResidue_ID_Finder
```
```bash
make 
```
- Or if you prefer manual building:

````bash
cmake -S . -B build
cmake --build build
```

This produces the executable `cactus_resid` inside the `build/` folder. You can copy it anywhere you like.

---

## Basic concepts

Before diving into usage it helps to know three terms:

**Residue** — one amino acid in the protein chain. Each residue has a three-letter code (`ASP`, `GLU`, `NAG`, etc.) and a sequence number that identifies its position in the chain.

**Carboxyl group** — the reactive tip of an aspartate or glutamate side chain, made of two oxygen atoms. In the PDB file these are named `OD1`/`OD2` (for ASP) and `OE1`/`OE2` (for GLU). This tool uses the midpoint between those two oxygens as the geometric probe for each residue.

**Cutoff distance** — a radius in Ångströms (Å). One Ångström is 0.1 nm; typical bond lengths are ~1.5 Å and typical interaction distances are 3–8 Å.

---

## Usage

```
cactus_resid <subcommand> [OPTIONS] [PATH]
```

`PATH` is the folder containing your PDB files. If you omit it the tool looks in the current directory.

There are two subcommands, described below. Run either one with `--help` for a full option list:

```bash
cactus_resid dist --help
cactus_resid coord --help
```

---

## Subcommand: `dist`

**What it does:** Finds pairs of acidic residues that are near a carbohydrate residue (such as a muramic acid sugar) in the same structure.

This is a two-stage search:

1. **Stage 1 — carbohydrate proximity.** The tool goes through every ASP and GLU in the structure and computes the distance from that residue's carboxyl midpoint to every atom of every carbohydrate residue in the file. If the midpoint is within `--carb-cutoff` Å of any carbohydrate atom, the residue becomes a *candidate*.

2. **Stage 2 — pair finding.** Among the candidates, every pair of distinct residues whose carboxyl groups are within `--pair-cutoff` Å of each other is reported.

### Options

| Flag | Default | Meaning |
|---|---|---|
| `--carb-cutoff DIST` | `8.0` | How close (Å) a carboxyl midpoint must be to a carbohydrate atom to qualify as a candidate |
| `--pair-cutoff DIST` | `20.0` | Maximum distance (Å) between two candidate residues to be reported as a pair |
| `--residues RES` | `GLU,ASP` | Which acidic residue types to search for |
| `--carbs RES` | *(built-in list)* | Extra carbohydrate residue names to add to the built-in list |
| `-o / --output FILE` | `results.log` | Where to write the results |

### Example

```bash
# Search all PDB files in ./structures/, looking for ASP/GLU pairs
# within 6 Å of any carbohydrate atom, where the pair members
# are within 15 Å of each other.
cactus_resid dist ./structures/ --carb-cutoff 6.0 --pair-cutoff 15.0

# Search for GLU pairs only, add a custom sugar code, save to a named file
cactus_resid dist ./structures/ --residues GLU --carbs MUR,NAM -o glu_pairs.log
```

### Built-in carbohydrate list

The tool recognises the following residue codes out of the box. You can extend this list at runtime with `--carbs`.

| Category | Codes |
|---|---|
| Hexoses | GLC, BGC, BGLC, GLA, GAL, MAN, BMA, FRU, TAG, SOR |
| Pentoses | RIB, ARA, XYL, LYX |
| Deoxy sugars | FUC, FCA, RHA, RHM |
| GlcNAc / GalNAc / ManNAc | NAG, NDG, NGA, A2G, BNG, MNA |
| Neuraminic acids | SIA, SLB, NGC |
| Disaccharides | LAT, MAL, CEL, TRE |
| Uronic acids | GCU, IDR |
| Other | SGN, SHB, GCS, LSB, LB |

---

## Subcommand: `coord`

**What it does:** The same two-stage search as `dist`, but instead of using a carbohydrate in the same file as the reference, it uses the atoms of a **separate reference PDB file** that you supply.

This is useful when you have a known structure — a co-crystallised ligand, an inhibitor, a catalytic residue from a well-characterised enzyme — and you want to find acidic residue pairs that sit in the geometrically equivalent position across a set of other structures.

### How it works

1. The reference PDB is parsed and a set of **anchor atoms** is extracted from it (you can filter by residue name and atom name — see below).
2. For each query PDB, the same two-stage pipeline runs with those anchor atoms as the reference set: any acidic residue whose carboxyl midpoint is within `--carb-cutoff` of any anchor atom becomes a candidate, and pairs within `--pair-cutoff` are reported.

### Options

| Flag | Default | Meaning |
|---|---|---|
| `--ref FILE` | *(required)* | Path to the reference PDB file |
| `--ref-residues RES` | *(all residues)* | Comma-separated residue names to extract from the reference |
| `--ref-atoms ATOMS` | *(all atoms)* | Comma-separated atom names to extract from those residues |
| `--ref-seqid IDS` | *(all sequence numbers)* | Comma-separated residue sequence numbers to extract from the reference |
| `--carb-cutoff DIST` | `8.0` | How close (Å) a carboxyl midpoint must be to a reference atom to qualify |
| `--pair-cutoff DIST` | `20.0` | Maximum distance (Å) between two candidates to form a pair |
| `--residues RES` | `GLU,ASP` | Which acidic residue types to search for in the query structures |
| `-o / --output FILE` | `coord_results.log` | Where to write the results |

### Examples

```bash
# Use every atom in ref.pdb as the anchor set
cactus_resid coord --ref ref.pdb ./structures/

# Use only the carboxyl oxygens of an ASP in the reference —
# this finds query residues that sit exactly where ASP's carboxyl sits
# in the reference structure
cactus_resid coord --ref ref.pdb \
    --ref-residues ASP --ref-atoms OD1,OD2 \
    ./structures/

# Use all atoms of a co-crystallised ligand (residue code LIG)
cactus_resid coord --ref ref.pdb --ref-residues LIG ./structures/

# Use carboxyl oxygens from both ASP and GLU in the reference
cactus_resid coord --ref ref.pdb \
    --ref-residues ASP,GLU \
    --ref-atoms OD1,OD2,OE1,OE2 \
    ./structures/

# Wider search, only look for GLU in the query structures
cactus_resid coord --ref ref.pdb \
    --ref-residues LIG \
    --carb-cutoff 12.0 --residues GLU \
    -o glu_near_lig.log \
    ./structures/

# Pin to one specific ASP by sequence number
cactus_resid coord --ref ref.pdb \
    --ref-residues ASP --ref-seqid 142 \
    ./structures/

# Use two specific residues by sequence number, only their carboxyl oxygens
cactus_resid coord --ref ref.pdb \
    --ref-residues ASP,GLU \
    --ref-seqid 142,158 \
    --ref-atoms OD1,OD2,OE1,OE2 \
    ./structures/
```

> **Note on multiple reference residues:** `--ref-residues` accepts a comma-separated list of any length. If your reference PDB contains several ASP residues and you write `--ref-residues ASP`, all of them are used as anchors. To pin the search to specific residues by position, use `--ref-seqid` with a comma-separated list of sequence numbers (e.g. `--ref-seqid 142,158`). The three reference filters — `--ref-residues`, `--ref-seqid`, and `--ref-atoms` — can be combined freely; an atom must pass every filter that is specified.

---

## Output format

Both subcommands write results to a log file (and print a progress bar to the terminal while running). The log is plain text and designed to be easy to search with `grep`.

### Log file structure

```
# [dist] Acidic-residue pairs near carbohydrate residues
# carb-cutoff : 8.0 Å
# pair-cutoff : 20.0 Å
# residues    : GLU, ASP
#
# FILE: ./structures/1abc.pdb
PAIR  file=./structures/1abc.pdb  res1=142  res2=158  dist=12.034
PAIR  file=./structures/1abc.pdb  res1=142  res2=201  dist=17.450
# FILE: ./structures/2xyz.pdb
  (no pairs found)
```

Lines beginning with `#` are comments describing the run parameters or the current file being processed. Lines beginning with `PAIR` are results.

### Recovering results with grep

```bash
# All pairs from any file
grep '^PAIR' results.log

# All pairs from a specific structure
grep '^PAIR.*1abc' results.log

# All pairs involving residue 142
grep 'res1=142\|res2=142' results.log

# Pairs closer than 10 Å (requires awk)
grep '^PAIR' results.log | awk -F'dist=' '{if ($2 < 10.0) print}'
```

---

## A note on distances

The tool does **not** apply periodic boundary conditions. This is intentional — it is designed for crystal structures and static models, not molecular dynamics trajectories. Applying periodic boundary conditions to crystal structures would produce spurious short-distance artefacts.

---

## Questions and common issues

**The tool reports no pairs for any file.**
Check that your PDB files contain HETATM records for carbohydrate residues (`dist`) or that your reference file contains the residue names you specified (`coord`). You can also try increasing `--carb-cutoff`.

**I have a sugar that isn't in the built-in list.**
Use `--carbs YOURCODE` to add it at runtime, e.g. `--carbs MUR`.

**I want to search only aspartates, not glutamates.**
Use `--residues ASP`.

**The progress bar looks garbled.**
This can happen if the terminal window is very narrow. It does not affect the results in the log file.

**I want to filter results by chain.**
The current output includes the file name and residue sequence numbers. Chain-level filtering is not yet supported in the output but can be done manually by cross-referencing the sequence numbers against the PDB file.
