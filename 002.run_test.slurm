#!/bin/bash
#SBATCH --job-name=DRef614_Equiv

set -euo pipefail

WORK_DIR="${SLURM_SUBMIT_DIR}"
source "${WORK_DIR}/000.config.sh"

echo "DOCK6.14 Dynamic Referencing Equivalence Test"
echo "System: $TEST_SYSTEM"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

validate_configuration

# Confirm that both calculations use identical input files
CANONICAL_HASH=$(sha256sum "$CANONICAL_INPUT" | awk '{print $1}')
BASELINE_HASH=$(sha256sum "$BASELINE_INPUT" | awk '{print $1}')
MODULARIZED_HASH=$(sha256sum "$MODULARIZED_INPUT" | awk '{print $1}')

if [[ "$CANONICAL_HASH" != "$BASELINE_HASH" ]] ||
   [[ "$CANONICAL_HASH" != "$MODULARIZED_HASH" ]]; then
    echo "ERROR: Input verification failed."
    exit 1
fi

# Run one DOCK6 implementation and record its completion status
run_dock()
{
    local LABEL="$1"
    local DOCK_BIN="$2"
    local RUN_DIR="$3"
    local INPUT_FILE="$4"
    local OUTPUT_FILE="$5"
    local LOG_FILE="$6"

    echo
    echo "Running $LABEL"
    echo "Executable: $DOCK_BIN"
    echo "Start: $(date)"

    (
        cd "$RUN_DIR"

        rm -f .success .failed

        if "$DOCK_BIN" -i "$INPUT_FILE" -o "$OUTPUT_FILE" > "$LOG_FILE" 2>&1; then
            if [[ -s "$OUTPUT_FILE" ]]; then
                touch .success
                echo "$LABEL: SUCCESS"
                exit 0
            fi
            echo "$LABEL: Missing DOCK output file"
        else
            echo "$LABEL: DOCK execution failed"
        fi

        touch .failed
        echo "$LABEL: FAILED"
        exit 1
    )
}

BASELINE_STATUS=0
MODULARIZED_STATUS=0

run_dock \
    "Baseline DOCK6.14" \
    "$BASELINE_DOCK_BIN" \
    "$BASELINE_RUN_DIR" \
    "$BASELINE_INPUT" \
    "$BASELINE_OUTPUT" \
    "$BASELINE_LOG" || BASELINE_STATUS=$?

run_dock \
    "Modularized DOCK6.14" \
    "$MODULARIZED_DOCK_BIN" \
    "$MODULARIZED_RUN_DIR" \
    "$MODULARIZED_INPUT" \
    "$MODULARIZED_OUTPUT" \
    "$MODULARIZED_LOG" || MODULARIZED_STATUS=$?

echo
echo "Baseline status: $BASELINE_STATUS"
echo "Modularized status: $MODULARIZED_STATUS"
echo "Finished: $(date)"

if (( BASELINE_STATUS != 0 || MODULARIZED_STATUS != 0 )); then
    exit 1
fi

echo "Both DOCK6.14 calculations completed."
