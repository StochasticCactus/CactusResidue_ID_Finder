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
git clone <git@github.com:StochasticCactus/CactusResidue_ID_Finder.git>
cd CactusResidue_ID_Finder
```
```bash
make 
```
- Or if you prefer manual building:

```bash
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

**What it does:** For each query PDB, finds the single acidic residue (ASP or GLU) whose carboxyl-group midpoint is nearest to each supplied target point in space. One result line is written per target point per file — either `MATCH` or `NOMATCH` — making it easy to recover residue numbers across large model sets (e.g. RFDiffusion outputs) where catalytic residues are geometrically fixed but their sequence numbers vary between models.

Only the carboxyl-group midpoint is used as the geometric probe:

| Residue | Oxygens used | Probe |
|---------|-------------|-------|
| ASP | OD1, OD2 | midpoint of OD1 and OD2 |
| GLU | OE1, OE2 | midpoint of OE1 and OE2 |

### Defining target points

You must supply at least one target point. Both sources may be combined freely; points are numbered in the order they are collected (reference residues first, then `--point` arguments).

**From a reference PDB** — provide a PDB file containing your known catalytic residues and list their sequence numbers. The tool computes the carboxyl midpoint of each and uses it as a target:

```bash
cactus_resid coord --ref ref.pdb --ref-seqid 142,158 ./models/
```

The reference residues must be ASP or GLU; other types are rejected with an error message.

**From explicit coordinates** — supply `--point x,y,z` directly on the command line. The flag is repeatable:

```bash
cactus_resid coord --point 5.0,2.0,7.0 --point 5.0,3.0,8.0 ./models/
```

**Mixed** — both sources at once:

```bash
cactus_resid coord --ref ref.pdb --ref-seqid 142 --point 5.0,3.0,8.0 ./models/
```

### Matching logic

For each target point the nearest acidic-residue midpoint in the query file is found. If its distance is within `--cutoff` the result is a `MATCH`; if no residue falls within the cutoff the result is `NOMATCH`. Exactly one line is written per target point per file — there is no ambiguity from multiple candidates.

### Options

| Flag | Default | Meaning |
|---|---|---|
| `--ref FILE` | — | Reference PDB from which to derive target points |
| `--ref-seqid IDS` | — | Comma-separated sequence numbers of the catalytic residues in the reference (required when `--ref` is given) |
| `--point x,y,z` | — | Explicit target coordinate in Å (repeatable) |
| `--cutoff DIST` | `5.0` | Max distance (Å) for a valid match; beyond this → `NOMATCH` |
| `--residues RES` | `GLU,ASP` | Acidic residue types to search for in the query structures |
| `-o / --output FILE` | `coord_results.log` | Where to write results |

### Output format

```
# [coord] Nearest acidic residue to each target point
# cutoff    : 5.0 Å
# residues  : GLU, ASP
# points    :
#   1  ref:ASP142  10.500  22.100  8.300
#   2  ref:GLU158  14.200  19.800  11.050
#
# FILE: ./models/design_0001.pdb
MATCH    file=./models/design_0001.pdb  point=1  res=87   resname=ASP  dist=1.203
MATCH    file=./models/design_0001.pdb  point=2  res=103  resname=GLU  dist=0.987
# FILE: ./models/design_0002.pdb
MATCH    file=./models/design_0002.pdb  point=1  res=91   resname=ASP  dist=2.441
NOMATCH  file=./models/design_0002.pdb  point=2  nearest=6.812  cutoff=5.0
```

`MATCH` lines report the residue sequence number, residue name, and distance. `NOMATCH` lines report how far away the nearest residue actually was, so you can judge whether to widen `--cutoff`.

### Filtering results with grep

```bash
# All successful matches
grep '^MATCH' coord_results.log

# All failures — models where a catalytic position was not recovered
grep '^NOMATCH' coord_results.log

# Just the residue numbers for point 1 across all models
grep '^MATCH.*point=1' coord_results.log | grep -oP 'res=\K[0-9]+'

# Files where point 2 failed
grep '^NOMATCH.*point=2' coord_results.log | grep -oP 'file=\K\S+'

# Models where BOTH points matched (no NOMATCH lines for that file)
awk '/^# FILE:/{f=$3} /^NOMATCH/{bad[f]=1} END{for(k in bad) print k}' \
    coord_results.log   # prints files with at least one NOMATCH
```

### Examples

```bash
# Derive both target points from a reference structure
cactus_resid coord --ref ref.pdb --ref-seqid 142,158 ./models/

# Explicit coordinates, tighter cutoff
cactus_resid coord \
    --point 10.5,22.1,8.3 \
    --point 14.2,19.8,11.05 \
    --cutoff 3.0 \
    ./models/

# One point from reference, one explicit; search GLU only
cactus_resid coord \
    --ref ref.pdb --ref-seqid 142 \
    --point 14.2,19.8,11.05 \
    --residues GLU \
    -o glu_search.log \
    ./models/
```

---

## Output format

Both subcommands write results to a log file (and print a progress bar to the terminal while running). The log is plain text and designed to be easy to search with `grep`.

### `dist` log file structure

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

### Recovering `dist` results with grep

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
