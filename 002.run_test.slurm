#!/bin/bash

#SBATCH --job-name=DRef_Equiv
#SBATCH --output=comparison_results/slurm_%j.out
#SBATCH --error=comparison_results/slurm_%j.err

set -euo pipefail

WORK_DIR="${SLURM_SUBMIT_DIR}"
source "${WORK_DIR}/000.config.sh"

cd "${WORK_DIR}"

echo
echo "DRef Equivalence Test"
echo
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "System: ${TEST_SYSTEM}"
echo "Anchor: ${TEST_ANCHOR}"
echo "Start: $(date)"
echo

validate_configuration

CANONICAL_HASH=$(sha256sum "${CANONICAL_INPUT}" | awk '{print $1}')
INHERITED_HASH=$(sha256sum "${INHERITED_INPUT}" | awk '{print $1}')
ISOLATED_HASH=$(sha256sum "${ISOLATED_INPUT}" | awk '{print $1}')

if [[ "${CANONICAL_HASH}" != "${INHERITED_HASH}" ]] ||
   [[ "${CANONICAL_HASH}" != "${ISOLATED_HASH}" ]]; then
    echo "ERROR: Equivalence-test inputs are not identical."
    exit 1
fi

run_dock()
{
    local LABEL="$1"
    local DOCK_BIN="$2"
    local RUN_DIR="$3"
    local INPUT_FILE="$4"
    local OUTPUT_FILE="$5"
    local LOG_FILE="$6"

    echo
    echo "Running ${LABEL}"
    echo "Executable: ${DOCK_BIN}"
    echo "Directory: ${RUN_DIR}"
    echo "Started: $(date)"
    echo

    cd "${RUN_DIR}"

    rm -f \
        .success \
        .failed \
        output*

    set +e

    "${DOCK_BIN}" \
        -i "${INPUT_FILE}" \
        -o "${OUTPUT_FILE}" \
        > "${LOG_FILE}" 2>&1

    DOCK_EXIT=$?

    set -e

    if [[ ${DOCK_EXIT} -eq 0 ]] &&
       [[ -s "${OUTPUT_FILE}" ]] &&
       grep -q "Total elapsed time" "${OUTPUT_FILE}"; then

        : > .success

        echo "${LABEL}: SUCCESS"
        echo "Finished: $(date)"

        return 0
    fi

    : > .failed

    echo "${LABEL}: FAILED"
    echo "Exit code: ${DOCK_EXIT}"
    echo "Finished: $(date)"

    return 1
}

INHERITED_STATUS=0
ISOLATED_STATUS=0

run_dock \
    "Inherited DRef" \
    "${INHERITED_DOCK_BIN}" \
    "${INHERITED_RUN_DIR}" \
    "${INHERITED_INPUT}" \
    "${INHERITED_OUTPUT}" \
    "${INHERITED_LOG}" || INHERITED_STATUS=$?

run_dock \
    "Isolated DRef" \
    "${ISOLATED_DOCK_BIN}" \
    "${ISOLATED_RUN_DIR}" \
    "${ISOLATED_INPUT}" \
    "${ISOLATED_OUTPUT}" \
    "${ISOLATED_LOG}" || ISOLATED_STATUS=$?

echo
echo "Execution Summary"
echo
echo "Inherited status: ${INHERITED_STATUS}"
echo "Isolated status:  ${ISOLATED_STATUS}"
echo "Node: $(hostname)"
echo "Finished: $(date)"
echo

if (( INHERITED_STATUS != 0 || ISOLATED_STATUS != 0 )); then
    echo "DRef equivalence calculation failed."
    exit 1
fi

echo "Both calculations completed successfully."
