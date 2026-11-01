#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/000.config.sh"

build_one() {
    local SRC="$1"
    local DEST="$2"
    local LABEL="$3"
    local MAKE_DIR

    if [[ -f "${SRC}/src/Makefile" ]]; then
        MAKE_DIR="${SRC}/src"
    elif [[ -f "${SRC}/Makefile" ]]; then
        MAKE_DIR="${SRC}"
    else
        echo "ERROR: Makefile not found for ${LABEL}" >&2
        return 1
    fi

    echo "Compiling ${LABEL} from ${SRC}"

    # Rebuild the selected DOCK6.14 implementation.
    make -C "${MAKE_DIR}" clean
    make -C "${MAKE_DIR}" -j 8

    if [[ ! -x "${SRC}/bin/dock6" ]]; then
        echo "ERROR: Executable not found: ${SRC}/bin/dock6" >&2
        return 1
    fi

    install -m 755 "${SRC}/bin/dock6" "${DEST}"

    echo "${LABEL}: ${DEST}"
}

[[ -d "${INHERITED_SRC}" ]] || {
    echo "ERROR: Inherited source directory missing." >&2
    exit 1
}

[[ -d "${ISOLATED_SRC}" ]] || {
    echo "ERROR: Isolated source directory missing." >&2
    exit 1
}

if [[ -e "${RUN_ROOT}" ]]; then
    echo "ERROR: Prepared experiment runs exist. Rebuilding would change binary provenance." >&2
    exit 1
fi

build_one "${INHERITED_SRC}" "${INHERITED_DOCK_BIN}" "Inherited DRef"
build_one "${ISOLATED_SRC}" "${ISOLATED_DOCK_BIN}" "Isolated DRef"

echo
echo "Compilation completed."

sha256sum "${INHERITED_DOCK_BIN}" "${ISOLATED_DOCK_BIN}"
