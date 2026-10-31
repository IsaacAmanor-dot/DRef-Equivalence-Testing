#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/000.config.sh"

INPUT_SOURCE="${1:-${CANONICAL_INPUT}}"

if [[ ! -f "$INPUT_SOURCE" ]]; then
    echo "ERROR: Input file not found: $INPUT_SOURCE"
    exit 1
fi

# Prevent overwriting previous calculations
for DIR in "$BASELINE_RUN_DIR" "$MODULARIZED_RUN_DIR"; do
    if [[ -d "$DIR" ]] && [[ -n "$(find "$DIR" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
        echo "ERROR: Run directory is not empty: $DIR"
        exit 1
    fi
done

# Construct the 1O86 input without changing the original template
python3 - "$INPUT_SOURCE" "$CANONICAL_INPUT" "$SPHERES" "$GRID_PREFIX" "$REF_LIGAND" << 'PY'
from pathlib import Path
import sys

source_path = Path(sys.argv[1])
destination = Path(sys.argv[2])

replacements = {
    "receptor_site_file": sys.argv[3],
    "descriptor_grid_score_grid_prefix": sys.argv[4],
    "descriptor_hms_score_ref_filename": sys.argv[5],
}

lines = source_path.read_text().splitlines(keepends=True)
counts = {key: 0 for key in replacements}
updated = []

for line in lines:
    fields = line.split()

    if fields and fields[0] in replacements:
        key = fields[0]
        counts[key] += 1
        newline = "\n" if line.endswith("\n") else ""
        updated.append(f"{key:<65} {replacements[key]}{newline}")
    else:
        updated.append(line)

if any(count != 1 for count in counts.values()):
    raise SystemExit(f"ERROR: System parameter counts are invalid: {counts}")

text = "".join(updated)

if "/002_DRef_Systems/1EVE/" in text or "/002_DRef_Systems/1EQH/" in text:
    raise SystemExit("ERROR: An old system-specific path remains in the input.")

destination.write_text(text)
print("Prepared 1O86 input:", destination)
PY

# Validate configuration before creating the run inputs
validate_configuration

mkdir -p "$BASELINE_RUN_DIR"
mkdir -p "$MODULARIZED_RUN_DIR"
mkdir -p "$RESULTS_ROOT"

cp "$CANONICAL_INPUT" "$BASELINE_INPUT"
cp "$CANONICAL_INPUT" "$MODULARIZED_INPUT"

CANONICAL_HASH=$(sha256sum "$CANONICAL_INPUT" | awk '{print $1}')
BASELINE_HASH=$(sha256sum "$BASELINE_INPUT" | awk '{print $1}')
MODULARIZED_HASH=$(sha256sum "$MODULARIZED_INPUT" | awk '{print $1}')

if [[ "$CANONICAL_HASH" != "$BASELINE_HASH" ]] ||
   [[ "$CANONICAL_HASH" != "$MODULARIZED_HASH" ]]; then
    echo "ERROR: Input files are not identical."
    exit 1
fi

echo
echo "DOCK6.14 DRef Equivalence Test"
echo "System: ${TEST_SYSTEM}"
echo "Reference: ${REF_LIGAND}"
echo "Spheres: ${SPHERES}"
echo "Grid: ${GRID_PREFIX}"
echo "Anchor: $(parameter_value dn_fraglib_anchor_file)"
echo
echo "Canonical SHA256:   ${CANONICAL_HASH}"
echo "Baseline SHA256:    ${BASELINE_HASH}"
echo "Modularized SHA256: ${MODULARIZED_HASH}"
echo
echo "Setup completed successfully."
