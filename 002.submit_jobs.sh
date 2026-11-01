#!/bin/bash
#SBATCH --job-name=DRef_equivalence
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/000.config.sh"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    [[ -s "${RUN_ROOT}/tasks.tsv" ]] || { echo 'ERROR: Run 001.setup_test.sh first' >&2; exit 1; }
    (( SLURM_CPUS > 0 && SLURM_CPUS <= 40 && MAX_CONCURRENT > 0 && MAX_CONCURRENT <= SLURM_CPUS )) || {
        echo 'ERROR: Invalid CPU or concurrency configuration' >&2; exit 1;
    }
    [[ -x "${INHERITED_DOCK_BIN}" && -x "${ISOLATED_DOCK_BIN}" ]] || {
        echo 'ERROR: Missing named DOCK executable' >&2; exit 1;
    }
    mkdir -p "${RESULTS_ROOT}"
    [[ ! -e "${RUN_ROOT}/submitted_jobs.tsv" ]] || {
        echo 'ERROR: Jobs already submitted; check submitted_jobs.tsv' >&2; exit 1;
    }
    cd "${SCRIPT_DIR}"
    : > "${RUN_ROOT}/submitted_jobs.tsv"
    for MODE in inherited isolated; do
        JOB_ID=$(sbatch --parsable \
            --partition="${SLURM_PARTITION}" \
            --nodes=1 --ntasks=1 --cpus-per-task="${SLURM_CPUS}" \
            --time="${SLURM_TIME}" --chdir="${SCRIPT_DIR}" \
            --job-name="DRef_${MODE}" \
            --output="${RESULTS_ROOT}/${MODE}_%j.out" \
            --error="${RESULTS_ROOT}/${MODE}_%j.err" \
            --export="ALL,DREF_MODE=${MODE}" "${SCRIPT_DIR}/002.submit_jobs.sh")
        printf '%s\t%s\n' "${MODE}" "${JOB_ID}" | tee -a "${RUN_ROOT}/submitted_jobs.tsv"
    done
    exit 0
fi

case "${DREF_MODE:-}" in
    inherited) DOCK_BIN="${INHERITED_DOCK_BIN}" ;;
    isolated) DOCK_BIN="${ISOLATED_DOCK_BIN}" ;;
    *) echo 'ERROR: Invalid DREF_MODE' >&2; exit 1 ;;
esac

[[ -x "$DOCK_BIN" ]] || { echo "ERROR: Missing executable: $DOCK_BIN" >&2; exit 1; }
[[ -s "${RUN_ROOT}/tasks.tsv" ]] || { echo 'ERROR: Missing task list' >&2; exit 1; }

# Protect against recompiling binaries after input preparation.
EXPECTED=$(awk -F '\t' -v mode="$DREF_MODE" '$1==mode {print $2}' "${RUN_ROOT}/binary_sha256.txt")
ACTUAL=$(sha256sum "$DOCK_BIN" | awk '{print $1}')
[[ -n "$EXPECTED" && "$EXPECTED" == "$ACTUAL" ]] || {
    echo 'ERROR: Executable hash changed after setup' >&2; exit 1;
}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
printf 'Mode: %s | Node: %s | Job: %s\n' "$DREF_MODE" "$(hostname)" "$SLURM_JOB_ID"

run_anchor() {
    local SYSTEM="$1" ANCHOR="$2" DIR RC=0
    DIR="${RUN_ROOT}/${DREF_MODE}/${SYSTEM}/anc_${ANCHOR}"
    [[ -s "${DIR}/dn_generic.in" ]] || { echo "ERROR: Missing input: $DIR" >&2; return 1; }
    [[ ! -e "${DIR}/.success" && ! -e "${DIR}/.failed" && ! -e "${DIR}/.running" ]] || {
        echo "ERROR: Prior run detected in $DIR" >&2; return 1;
    }
    (
        cd "$DIR"
        touch .running
        "$DOCK_BIN" -i dn_generic.in -o dn_generic.out > dock_stdout.log 2>&1 || RC=$?
        rm -f .running
        if (( RC == 0 )) && [[ -s dn_generic.out ]]; then
            touch .success
            echo "Completed ${DREF_MODE}: ${SYSTEM} anc_${ANCHOR}"
        else
            touch .failed
            echo "FAILED ${DREF_MODE}: ${SYSTEM} anc_${ANCHOR}, exit=${RC}" >&2
            exit 1
        fi
    )
}

PIDS=()
FAILED=0
while IFS=$'\t' read -r SYSTEM ANCHOR; do
    [[ -n "$SYSTEM" ]] || continue
    run_anchor "$SYSTEM" "$ANCHOR" &
    PIDS+=("$!")
    if (( ${#PIDS[@]} >= MAX_CONCURRENT )); then
        for PID in "${PIDS[@]}"; do wait "$PID" || FAILED=1; done
        PIDS=()
    fi
done < "${RUN_ROOT}/tasks.tsv"
for PID in "${PIDS[@]}"; do wait "$PID" || FAILED=1; done
(( FAILED == 0 )) || { echo "ERROR: Some ${DREF_MODE} calculations failed" >&2; exit 1; }
echo "All ${DREF_MODE} calculations completed"
