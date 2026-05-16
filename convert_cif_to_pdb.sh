#!/usr/bin/env bash
# =============================================================================
# convert_cif_to_pdb.sh
#
# Recursively converts all .cif.gz files to .pdb using Open Babel.
# Mirrors the subdirectory structure into output_dir.
# macOS-compatible (bash 3.2, no flock, no mapfile, no Linux mktemp flags).
#
# Usage:
#   chmod +x convert_cif_to_pdb.sh
#   ./convert_cif_to_pdb.sh [OPTIONS] <input_dir> [output_dir]
#
# Options:
#   -j N   Parallel jobs (default: all CPU cores)
#   -f     Force overwrite existing .pdb files
#   -v     Verbose: print obabel command for each file
#   -h     Show this help
#
# Examples:
#   ./convert_cif_to_pdb.sh output_RFD3
#   ./convert_cif_to_pdb.sh -j 4 -f output_RFD3 ./pdbs_RFD3
# =============================================================================

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
JOBS=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
FORCE=false
VERBOSE=false
INPUT_DIR=""
OUTPUT_DIR=""

# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------
if [ -t 1 ]; then
    RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
    CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'
else
    RED=''; GREEN=''; YELLOW=''; CYAN=''; BOLD=''; RESET=''
fi

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
err()     { echo -e "${RED}[ERROR]${RESET} $*"; }

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------
while getopts ":j:fvh" opt; do
    case $opt in
        j) JOBS="$OPTARG" ;;
        f) FORCE=true ;;
        v) VERBOSE=true ;;
        h) sed -n '3,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        :) echo "Option -$OPTARG requires an argument."; exit 1 ;;
       \?) echo "Unknown option: -$OPTARG"; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

INPUT_DIR="${1:-}"
OUTPUT_DIR="${2:-}"

if [[ -z "$INPUT_DIR" ]]; then
    echo "Usage: $0 [OPTIONS] <input_dir> [output_dir]"; exit 1
fi
if [[ ! -d "$INPUT_DIR" ]]; then
    err "Not a directory: $INPUT_DIR"; exit 1
fi

INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"

if [[ -n "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
    OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
fi

# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------
if ! command -v obabel &>/dev/null; then
    err "obabel not found. Install with:  brew install open-babel"
    exit 1
fi
info "obabel: $(obabel --version 2>&1 | grep -o 'Open Babel [0-9.]*' || echo 'found')"

# ---------------------------------------------------------------------------
# Temp workspace (macOS-safe: no suffix flags, no flock)
# ---------------------------------------------------------------------------
WORK_DIR=$(mktemp -d)
LOG_DIR="$WORK_DIR/logs"
TMP_CIF_DIR="$WORK_DIR/cifs"
mkdir -p "$LOG_DIR" "$TMP_CIF_DIR"

cleanup() { rm -rf "$WORK_DIR"; }
trap cleanup EXIT

# ---------------------------------------------------------------------------
# Collect files (NUL-delimited, bash 3.2 safe — no mapfile)
# ---------------------------------------------------------------------------
FILELIST="$WORK_DIR/filelist"
find "$INPUT_DIR" -type f -iname "*.cif.gz" -print0 > "$FILELIST"
TOTAL=$(tr -cd '\0' < "$FILELIST" | wc -c | tr -d ' ')

if [[ "$TOTAL" -eq 0 ]]; then
    warn "No .cif.gz files found under: $INPUT_DIR"; exit 0
fi

info "Found ${BOLD}${TOTAL}${RESET} .cif.gz file(s)"
info "Parallel jobs  : $JOBS"
info "Force overwrite: $FORCE"
[[ -n "$OUTPUT_DIR" ]] \
    && info "Output dir     : $OUTPUT_DIR" \
    || info "Output dir     : in-place (next to each .cif.gz)"
echo ""

# ---------------------------------------------------------------------------
# Per-file worker — called by xargs
# Each job writes a one-line result to its own log file (no locking needed)
# ---------------------------------------------------------------------------
convert_one() {
    local cif_gz="$1"
    local input_dir="$2"
    local output_dir="$3"
    local force="$4"
    local verbose="$5"
    local log_dir="$6"
    local tmp_cif_dir="$7"

    local base src_dir dest_dir rel_dir pdb_out tmp_cif log_file

    base=$(basename "$cif_gz" .cif.gz)
    src_dir=$(dirname "$cif_gz")

    # Mirror subdirectory structure
    if [[ -n "$output_dir" ]]; then
        rel_dir="${src_dir#"$input_dir"}"
        rel_dir="${rel_dir#/}"
        dest_dir="${output_dir}${rel_dir:+/$rel_dir}"
    else
        dest_dir="$src_dir"
    fi

    pdb_out="$dest_dir/${base}.pdb"

    # Unique log file per job (PID + RANDOM — no flock needed)
    log_file="$log_dir/$$.${RANDOM}.log"

    # Skip?
    if [[ -f "$pdb_out" && "$force" == "false" ]]; then
        echo "SKIP" > "$log_file"
        echo "[SKIP]  $pdb_out"
        return
    fi

    mkdir -p "$dest_dir"

    # Unique temp .cif — PID + RANDOM avoids collisions; no suffix needed by mktemp
    tmp_cif="${tmp_cif_dir}/${$}_${RANDOM}"

    if ! gunzip -c "$cif_gz" > "$tmp_cif" 2>/dev/null; then
        echo "FAIL" > "$log_file"
        err "Decompress failed: $cif_gz"
        rm -f "$tmp_cif"
        return
    fi

    if [[ "$verbose" == "true" ]]; then
        echo "[CMD]   obabel -i cif $tmp_cif -o pdb -O $pdb_out"
    fi

    if obabel -i cif "$tmp_cif" -o pdb -O "$pdb_out" 2>/dev/null; then
        if [[ -s "$pdb_out" ]]; then
            echo "OK" > "$log_file"
            echo "[OK]    $pdb_out"
        else
            echo "FAIL" > "$log_file"
            err "Empty output for: $cif_gz"
            rm -f "$pdb_out"
        fi
    else
        echo "FAIL" > "$log_file"
        err "obabel failed for: $cif_gz"
        rm -f "$pdb_out"
    fi

    rm -f "$tmp_cif"
}

export -f convert_one

# ---------------------------------------------------------------------------
# Run in parallel
# ---------------------------------------------------------------------------
xargs -0 -P "$JOBS" -I {} bash -c \
    'convert_one "$@"' _ \
    {} \
    "$INPUT_DIR" \
    "$OUTPUT_DIR" \
    "$FORCE" \
    "$VERBOSE" \
    "$LOG_DIR" \
    "$TMP_CIF_DIR" \
    < "$FILELIST"

# ---------------------------------------------------------------------------
# Tally results from log files (no shared counters, no flock needed)
# ---------------------------------------------------------------------------
OK=$(grep -rl   '^OK$'   "$LOG_DIR" 2>/dev/null | wc -l | tr -d ' ')
SKIP=$(grep -rl '^SKIP$' "$LOG_DIR" 2>/dev/null | wc -l | tr -d ' ')
FAIL=$(grep -rl '^FAIL$' "$LOG_DIR" 2>/dev/null | wc -l | tr -d ' ')

echo ""
echo -e "${BOLD}=============================${RESET}"
echo -e "${BOLD} Conversion summary${RESET}"
echo -e "${BOLD}=============================${RESET}"
echo -e "  Total        : $TOTAL"
echo -e "  ${GREEN}Converted OK ${RESET}: $OK"
echo -e "  ${YELLOW}Skipped      ${RESET}: $SKIP"
echo -e "  ${RED}Failed       ${RESET}: $FAIL"
echo ""

if [[ "$FAIL" -gt 0 ]]; then
    warn "$FAIL file(s) failed."
    exit 1
fi

success "All done."
