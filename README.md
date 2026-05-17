## Subcommand: `coord`

**What it does:** For each query PDB, finds the single acidic residue (ASP or
GLU) whose carboxyl-group midpoint is nearest to each supplied target point in
space. One result line is written per target point per file — either `MATCH` or
`NOMATCH` — making it easy to recover residue numbers across large model sets
(e.g. RFDiffusion outputs) where catalytic residues are geometrically fixed but
their sequence numbers vary between models.

Only the carboxyl-group midpoint is used as the geometric probe:

| Residue | Oxygens used | Probe |
|---------|-------------|-------|
| ASP | OD1, OD2 | midpoint of OD1 and OD2 |
| GLU | OE1, OE2 | midpoint of OE1 and OE2 |

### Defining target points

You must supply at least one target point. Both sources may be combined freely;
points are numbered in the order they are collected (reference residues first,
then `--point` arguments).

**From a reference PDB** — provide a PDB file containing your known catalytic
residues and list their sequence numbers. The tool computes the carboxyl
midpoint of each and uses it as a target:

```bash
cactus_resid coord --ref ref.pdb --ref-seqid 142,158 ./models/
```

The reference residues must be ASP or GLU; other types are rejected with an
error message.

**From explicit coordinates** — supply `--point x,y,z` directly on the command
line. The flag is repeatable:

```bash
cactus_resid coord --point 5.0,2.0,7.0 --point 5.0,3.0,8.0 ./models/
```

**Mixed** — both sources at once:

```bash
cactus_resid coord --ref ref.pdb --ref-seqid 142 --point 5.0,3.0,8.0 ./models/
```

### Matching logic

For each target point the nearest acidic-residue midpoint in the query file is
found. If its distance is within `--cutoff` the result is a `MATCH`; if no
residue falls within the cutoff the result is `NOMATCH`. Exactly one line is
written per target point per file — there is no ambiguity from multiple
candidates.

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

`MATCH` lines report the residue sequence number, residue name, and distance.
`NOMATCH` lines report how far away the nearest residue actually was, so you
can judge whether to widen `--cutoff`.

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
