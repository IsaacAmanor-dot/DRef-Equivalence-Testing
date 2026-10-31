#!/bin/bash
set -euo pipefail

ROOT="/gpfs/projects/rizzo/iamanor/Collaborations/Carissa/002_Dref_Isolated_6.14_version"
MOD_SOURCE="${ROOT}/dock6_dev"
BASE_SOURCE="${ROOT}/dock6_baseline_source_7d38760d"
BUILD_DIR="${ROOT}/gpu_long_builds"
BASE_COMMIT="7d38760dacb1b4a1c8b28e6eaaeab4a833fc89a6"

# Verify the compute-node allocation
if [[ "${SLURM_JOB_PARTITION:-}" != "gpu-long" ]]; then
    echo "ERROR: Run inside a gpu-long allocation."
    exit 1
fi

module load gcc/13.2.0
module load git

mkdir -p "$BUILD_DIR"

# Verify the modularized source
if [[ "$(git -C "$MOD_SOURCE" branch --show-current)" != "dref-modularization" ]]; then
    echo "ERROR: Expected dref-modularization branch."
    exit 1
fi

if [[ ! -f "${MOD_SOURCE}/src/dock/dynamic_reference.cpp" ||
      ! -f "${MOD_SOURCE}/src/dock/dynamic_reference.h" ]]; then
    echo "ERROR: Dynamic Reference source files are missing."
    exit 1
fi

# Prepare the original DOCK6.14 source at the baseline commit
if [[ ! -e "$BASE_SOURCE" ]]; then
    git -C "$MOD_SOURCE" worktree add --detach "$BASE_SOURCE" "$BASE_COMMIT"
fi

if [[ "$(git -C "$BASE_SOURCE" rev-parse HEAD)" != "$BASE_COMMIT" ]]; then
    echo "ERROR: Incorrect baseline source commit."
    exit 1
fi

if ! git -C "$BASE_SOURCE" diff --quiet ||
   ! git -C "$BASE_SOURCE" diff --cached --quiet; then
    echo "ERROR: Baseline tracked source has modifications."
    exit 1
fi

# Record the compilation environment
hostname > "${BUILD_DIR}/build_node.txt"
echo "$SLURM_JOB_PARTITION" > "${BUILD_DIR}/build_partition.txt"
g++ --version | head -1 > "${BUILD_DIR}/compiler_version.txt"

git -C "$MOD_SOURCE" rev-parse HEAD > "${BUILD_DIR}/modularized_source_commit.txt"
git -C "$MOD_SOURCE" status --short > "${BUILD_DIR}/modularized_source_status.txt"

# Compile an implementation and preserve its executable
build_dock()
{
    local LABEL="$1"
    local SOURCE="$2"
    local OUTPUT="$3"

    echo "Compiling $LABEL on $(hostname)"

    cd "${SOURCE}/install"

    ./configure gnu > "${BUILD_DIR}/${LABEL}_configure.log" 2>&1
    make cleandock > "${BUILD_DIR}/${LABEL}_clean.log" 2>&1

    make dock -j 8 > "${BUILD_DIR}/${LABEL}_build.log" 2>&1 || {
        tail -80 "${BUILD_DIR}/${LABEL}_build.log"
        exit 1
    }

    cd "$SOURCE"

    if [[ ! -x bin/dock6 ]]; then
        echo "ERROR: $LABEL executable was not produced."
        exit 1
    fi

    cp -p bin/dock6 "$OUTPUT"
    sha256sum "$OUTPUT" > "${BUILD_DIR}/${LABEL}.sha256"

    echo "$LABEL compilation successful."
}

build_dock "baseline" "$BASE_SOURCE" "${BUILD_DIR}/dock6.baseline.gpu"
build_dock "modularized" "$MOD_SOURCE" "${BUILD_DIR}/dock6.modularized.gpu"

ls -lh "${BUILD_DIR}/dock6.baseline.gpu" "${BUILD_DIR}/dock6.modularized.gpu"

cat "${BUILD_DIR}/baseline.sha256"
cat "${BUILD_DIR}/modularized.sha256"

echo "Both DOCK6.14 implementations compiled successfully."
