#!/bin/bash

#SBATCH --job-name=DRef_equivalence

set -euo pipefail

if [[ "${1:-submit}" == "submit" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
else
    SCRIPT_DIR="${SLURM_SUBMIT_DIR:?Missing SLURM_SUBMIT_DIR}"
fi
source "${SCRIPT_DIR}/000.config.sh"

MODE="${1:-submit}"

if [[ "$MODE" == "submit" ]]; then

    [[ -s "${RUN_ROOT}/tasks.tsv" ]] || {
        echo "ERROR: Run 001.setup_test.sh first" >&2
        exit 1
    }

    (( SLURM_CPUS > 0 && SLURM_CPUS <= 40 )) || {
        echo "ERROR: Invalid CPU allocation" >&2
        exit 1
    }

    (( MAX_CONCURRENT > 0 && MAX_CONCURRENT <= SLURM_CPUS )) || {
        echo "ERROR: Invalid concurrency setting" >&2
        exit 1
    }

    [[ -x "${INHERITED_DOCK_BIN}" ]] || {
        echo "ERROR: Inherited executable missing" >&2
        exit 1
    }

    [[ -x "${ISOLATED_DOCK_BIN}" ]] || {
        echo "ERROR: Isolated executable missing" >&2
        exit 1
    }

    [[ -s "${RUN_ROOT}/binary_sha256.txt" ]] || {
        echo "ERROR: Binary checksum file missing" >&2
        exit 1
    }

    [[ ! -e "${RUN_ROOT}/submitted_jobs.tsv" ]] || {
        echo "ERROR: Jobs previously submitted" >&2
        exit 1
    }

    mkdir -p "${RESULTS_ROOT}"

    for JOB_MODE in inherited isolated; do

        JOB_ID=$(sbatch --parsable \
            --partition="${SLURM_PARTITION}" \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task="${SLURM_CPUS}" \
            --time="${SLURM_TIME}" \
            --chdir="${SCRIPT_DIR}" \
            --job-name="DRef_${JOB_MODE}" \
            --output="${RESULTS_ROOT}/${JOB_MODE}_%j.out" \
            --error="${RESULTS_ROOT}/${JOB_MODE}_%j.err" \
            --export=ALL \
            "${SCRIPT_DIR}/002.submit_jobs.sh" "${JOB_MODE}")

        printf '%s\t%s\n' "${JOB_MODE}" "${JOB_ID}" \
            | tee -a "${RUN_ROOT}/submitted_jobs.tsv"

    done

    exit 0
fi

case "$MODE" in
    inherited)
        DOCK_BIN="${INHERITED_DOCK_BIN}"
        ;;
    isolated)
        DOCK_BIN="${ISOLATED_DOCK_BIN}"
        ;;
    *)
        echo "ERROR: Invalid experiment mode: ${MODE}" >&2
        exit 1
        ;;
esac

[[ -n "${SLURM_JOB_ID:-}" ]] || {
    echo "ERROR: Execution requires a SLURM allocation" >&2
    exit 1
}

[[ -x "$DOCK_BIN" ]] || {
    echo "ERROR: Missing executable: $DOCK_BIN" >&2
    exit 1
}

[[ -s "${RUN_ROOT}/tasks.tsv" ]] || {
    echo "ERROR: Missing task list" >&2
    exit 1
}

# Verify that the executable has not changed since setup.
EXPECTED=$(awk -F '\t' -v mode="$MODE" '$1 == mode {print $2}' "${RUN_ROOT}/binary_sha256.txt")
ACTUAL=$(sha256sum "$DOCK_BIN" | awk '{print $1}')

[[ -n "$EXPECTED" && "$EXPECTED" == "$ACTUAL" ]] || {
    echo "ERROR: Executable checksum mismatch for ${MODE}" >&2
    exit 1
}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Mode: ${MODE}"
echo "Executable: ${DOCK_BIN}"
echo "Node: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "Started: $(date)"

run_anchor() {

    local SYSTEM="$1"
    local ANCHOR="$2"
    local DIR="${RUN_ROOT}/${MODE}/${SYSTEM}/anc_${ANCHOR}"

    local RC=0

    [[ -s "${DIR}/dn_generic.in" ]] || {
        echo "ERROR: Missing input: ${DIR}" >&2
        return 1
    }

    [[ ! -e "${DIR}/.success" &&
       ! -e "${DIR}/.failed" &&
       ! -e "${DIR}/.running" ]] || {
        echo "ERROR: Previous calculation exists: ${DIR}" >&2
        return 1
    }

    (
        cd "${DIR}"

        touch .running

        # Execute one anchor independently.
        "${DOCK_BIN}" -i dn_generic.in -o dn_generic.out \
            > dock_stdout.log 2>&1 || RC=$?

        rm -f .running

        if (( RC == 0 )) && [[ -s dn_generic.out ]]; then
            touch .success
            echo "Completed ${MODE}: ${SYSTEM} anc_${ANCHOR}"
        else
            touch .failed
            echo "FAILED ${MODE}: ${SYSTEM} anc_${ANCHOR}" >&2
            exit 1
        fi
    )
}

PIDS=()
FAILED=0

# Process anchors concurrently within the allocated node.
while IFS=$'\t' read -r SYSTEM ANCHOR; do

    [[ -n "$SYSTEM" ]] || continue

    run_anchor "$SYSTEM" "$ANCHOR" &
    PIDS+=("$!")

    if (( ${#PIDS[@]} >= MAX_CONCURRENT )); then
        for PID in "${PIDS[@]}"; do
            wait "$PID" || FAILED=1
        done
        PIDS=()
    fi

done < "${RUN_ROOT}/tasks.tsv"

for PID in "${PIDS[@]}"; do
    wait "$PID" || FAILED=1
done

echo "Finished: $(date)"

if (( FAILED != 0 )); then
    echo "ERROR: Some ${MODE} calculations failed" >&2
    exit 1
fi

echo "All ${MODE} calculations completed"
