#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/000.config.sh"

# Validate the DOCK6.14 executables and system inputs
validate_configuration

# Verify the prepared input files
for FILE in "$CANONICAL_INPUT" "$BASELINE_INPUT" "$MODULARIZED_INPUT"; do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing input file: $FILE"
        exit 1
    fi
done

CANONICAL_HASH=$(sha256sum "$CANONICAL_INPUT" | awk '{print $1}')
BASELINE_HASH=$(sha256sum "$BASELINE_INPUT" | awk '{print $1}')
MODULARIZED_HASH=$(sha256sum "$MODULARIZED_INPUT" | awk '{print $1}')

if [[ "$CANONICAL_HASH" != "$BASELINE_HASH" ]] ||
   [[ "$CANONICAL_HASH" != "$MODULARIZED_HASH" ]]; then
    echo "ERROR: Input files are not identical."
    exit 1
fi

# Preserve completed calculations
for DIR in "$BASELINE_RUN_DIR" "$MODULARIZED_RUN_DIR"; do
    if [[ -e "$DIR/.success" || -e "$DIR/.failed" ]]; then
        echo "ERROR: Previous calculation markers exist in $DIR"
        exit 1
    fi
done

mkdir -p "$RESULTS_ROOT"

# Submit the baseline and modularized calculations on one node
JOB_ID=$(sbatch \
    --parsable \
    --chdir="$SCRIPT_DIR" \
    --partition="$SLURM_PARTITION" \
    --nodes="$SLURM_NODES" \
    --ntasks="$SLURM_NTASKS" \
    --time="$SLURM_TIME" \
    --output="$RESULTS_ROOT/slurm_%j.out" \
    --error="$RESULTS_ROOT/slurm_%j.err" \
    "$SCRIPT_DIR/002.run_test.slurm")

echo "DOCK6.14 DRef equivalence test submitted."
echo "System: $TEST_SYSTEM"
echo "Job ID: $JOB_ID"
echo "Input SHA256: $CANONICAL_HASH"
echo
echo "Check job: squeue -j $JOB_ID"
