#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/000.config.sh"

for FILE in \
    "${CANONICAL_INPUT}" \
    "${INHERITED_INPUT}" \
    "${ISOLATED_INPUT}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Missing input: ${FILE}"
        exit 1
    fi
done

CANONICAL_HASH=$(sha256sum "${CANONICAL_INPUT}" | awk '{print $1}')
INHERITED_HASH=$(sha256sum "${INHERITED_INPUT}" | awk '{print $1}')
ISOLATED_HASH=$(sha256sum "${ISOLATED_INPUT}" | awk '{print $1}')

echo
echo "Input verification"
echo
echo "Canonical: ${CANONICAL_HASH}"
echo "Inherited: ${INHERITED_HASH}"
echo "Isolated:  ${ISOLATED_HASH}"
echo

if [[ "${CANONICAL_HASH}" != "${INHERITED_HASH}" ]] ||
   [[ "${CANONICAL_HASH}" != "${ISOLATED_HASH}" ]]; then

    echo "ERROR: Equivalence-test inputs are not identical."
    echo "Run 001.setup_test.sh before submission."
    exit 1
fi

JOB_ID=$(sbatch \
    --parsable \
    --partition="${SLURM_PARTITION}" \
    --nodes="${SLURM_NODES}" \
    --ntasks="${SLURM_NTASKS}" \
    --time="${SLURM_TIME}" \
    "${SCRIPT_DIR}/002.run_test.slurm")

echo "DRef equivalence test submitted."
echo "Job ID: ${JOB_ID}"
echo
echo "Check status:"
echo "squeue -j ${JOB_ID}"
echo
