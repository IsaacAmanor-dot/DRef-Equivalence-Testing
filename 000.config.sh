#!/bin/bash

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_ROOT="${CONFIG_DIR}"

# Test system and anchor

TEST_SYSTEM="1EQH"

SYSTEM_ROOT="/gpfs/projects/rizzo/iamanor/Systems_and_Library_Files/001_Systems_Files/002_DRef_Systems"

ANCHOR_ROOT="/gpfs/projects/rizzo/iamanor/DOCK6_Development/Dynamics_Referencing/for_isaac/Dynamic_Reference_for_Isaac/anchors"

TEST_ANCHOR="anchor_1.mol2"
ANCHOR_FILE="${ANCHOR_ROOT}/${TEST_ANCHOR}"

# DOCK6.9 implementations

INHERITED_DOCK_ROOT="/gpfs/projects/rizzo/iamanor/DOCK6_Development/Dynamics_Referencing/for_isaac/Dynamic_Reference_for_Isaac/dock6_dref69_inherited"

INHERITED_DOCK_BIN="${INHERITED_DOCK_ROOT}/bin/dock6.dref_inherited"

ISOLATED_DOCK_ROOT="/gpfs/projects/rizzo/iamanor/DOCK6_Development/Dynamics_Referencing/for_isaac/Dynamic_Reference_for_Isaac/dock6.21_11_05_Dynamic_Reference_Developmental_V1.6.2"

ISOLATED_DOCK_BIN="${ISOLATED_DOCK_ROOT}/bin/dock6.dref_isolated"

# DOCK parameter files

DOCK_PARAMS="${ISOLATED_DOCK_ROOT}/parameters"

VDW_DEFN_FILE="${DOCK_PARAMS}/vdw_de_novo.defn"

FLEX_DEFN_FILE="/gpfs/projects/rizzo/zzz.programs/dock6.9_mpiv2018.0.3/parameters/flex.defn"

FLEX_DRIVE_FILE="/gpfs/projects/rizzo/zzz.programs/dock6.9_mpiv2018.0.3/parameters/flex_drive.tbl"

CHEM_DEFN_FILE="${DOCK_PARAMS}/chem.defn"

# De Novo fragment library

FRAGLIB_ROOT="/gpfs/projects/rizzo/iamanor/Systems_and_Library_Files/003_Libraries/DOCK_DN_Generic_Library/DOCK6.13_Library"

FRAGLIB_SCAFFOLD="${FRAGLIB_ROOT}/fraglib_scaffold.mol2"
FRAGLIB_LINKER="${FRAGLIB_ROOT}/fraglib_linker.mol2"
FRAGLIB_SIDECHAIN="${FRAGLIB_ROOT}/fraglib_sidechain.mol2"
FRAGLIB_TORENV="${FRAGLIB_ROOT}/fraglib_torenv.dat"

# Receptor-system files

SYSTEM_DIR="${SYSTEM_ROOT}/${TEST_SYSTEM}"

REF_LIGAND="${SYSTEM_DIR}/${TEST_SYSTEM}.lig.am1bcc.mol2"
SPHERES="${SYSTEM_DIR}/${TEST_SYSTEM}.rec.clust.close.sph"

GRID_PREFIX="${SYSTEM_DIR}/${TEST_SYSTEM}.rec"
GRID_NRG="${GRID_PREFIX}.nrg"
GRID_BMP="${GRID_PREFIX}.bmp"

# Runtime directories

INHERITED_ROOT="${WORK_ROOT}/inherited"
ISOLATED_ROOT="${WORK_ROOT}/isolated"

INHERITED_RUN_DIR="${INHERITED_ROOT}/system"
ISOLATED_RUN_DIR="${ISOLATED_ROOT}/system"

RESULTS_ROOT="${WORK_ROOT}/comparison_results"

mkdir -p \
    "${INHERITED_RUN_DIR}" \
    "${ISOLATED_RUN_DIR}" \
    "${RESULTS_ROOT}"

# Input and output files

CANONICAL_INPUT="${WORK_ROOT}/dn_generic.in"

INHERITED_INPUT="${INHERITED_RUN_DIR}/dn_generic.in"
ISOLATED_INPUT="${ISOLATED_RUN_DIR}/dn_generic.in"

INHERITED_OUTPUT="${INHERITED_RUN_DIR}/dn_generic.out"
ISOLATED_OUTPUT="${ISOLATED_RUN_DIR}/dn_generic.out"

INHERITED_LOG="${RESULTS_ROOT}/inherited_run.log"
ISOLATED_LOG="${RESULTS_ROOT}/isolated_run.log"

COMPARISON_REPORT="${RESULTS_ROOT}/equivalence_report.txt"

# SLURM resources

SLURM_PARTITION="rn-long-40core"
SLURM_TIME="2-00:00:00"
SLURM_NODES=1
SLURM_NTASKS=40

check_file()
{
    local NAME="$1"
    local FILE="$2"

    if [[ ! -f "${FILE}" ]]; then
        echo "ERROR: Missing ${NAME}: ${FILE}"
        return 1
    fi
}

check_executable()
{
    local NAME="$1"
    local FILE="$2"

    if [[ ! -x "${FILE}" ]]; then
        echo "ERROR: Missing ${NAME}: ${FILE}"
        return 1
    fi
}

validate_configuration()
{
    local ERRORS=0

    check_executable "inherited DOCK executable" "${INHERITED_DOCK_BIN}" \
        || ERRORS=$((ERRORS + 1))

    check_executable "isolated DOCK executable" "${ISOLATED_DOCK_BIN}" \
        || ERRORS=$((ERRORS + 1))

    check_file "reference ligand" "${REF_LIGAND}" \
        || ERRORS=$((ERRORS + 1))

    check_file "receptor spheres" "${SPHERES}" \
        || ERRORS=$((ERRORS + 1))

    check_file "grid NRG file" "${GRID_NRG}" \
        || ERRORS=$((ERRORS + 1))

    check_file "grid BMP file" "${GRID_BMP}" \
        || ERRORS=$((ERRORS + 1))

    check_file "anchor" "${ANCHOR_FILE}" \
        || ERRORS=$((ERRORS + 1))

    check_file "scaffold library" "${FRAGLIB_SCAFFOLD}" \
        || ERRORS=$((ERRORS + 1))

    check_file "linker library" "${FRAGLIB_LINKER}" \
        || ERRORS=$((ERRORS + 1))

    check_file "sidechain library" "${FRAGLIB_SIDECHAIN}" \
        || ERRORS=$((ERRORS + 1))

    check_file "torsion environment table" "${FRAGLIB_TORENV}" \
        || ERRORS=$((ERRORS + 1))

    check_file "VDW definition" "${VDW_DEFN_FILE}" \
        || ERRORS=$((ERRORS + 1))

    check_file "flex definition" "${FLEX_DEFN_FILE}" \
        || ERRORS=$((ERRORS + 1))

    check_file "flex drive table" "${FLEX_DRIVE_FILE}" \
        || ERRORS=$((ERRORS + 1))

    if (( ERRORS > 0 )); then
        echo
        echo "Configuration validation failed with ${ERRORS} error(s)."
        return 1
    fi

    return 0
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then

    echo
    echo "DRef Equivalence Testing"
    echo
    echo "System:              ${TEST_SYSTEM}"
    echo "Anchor:              ${ANCHOR_FILE}"
    echo "Inherited DOCK:      ${INHERITED_DOCK_BIN}"
    echo "Isolated DOCK:       ${ISOLATED_DOCK_BIN}"
    echo "Reference ligand:    ${REF_LIGAND}"
    echo "Spheres:             ${SPHERES}"
    echo "Grid prefix:         ${GRID_PREFIX}"
    echo "Fragment library:    ${FRAGLIB_ROOT}"
    echo "Partition:           ${SLURM_PARTITION}"
    echo

    validate_configuration

    echo
    echo "Configuration validation passed."
    echo

fi
