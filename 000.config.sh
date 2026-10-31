#!/bin/bash

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_ROOT="${CONFIG_DIR}"

# Test system
TEST_SYSTEM="1O86"
SYSTEM_ROOT="/gpfs/projects/rizzo/iamanor/Systems_and_Library_Files/001_Systems_Files/002_DRef_Systems"
SYSTEM_DIR="${SYSTEM_ROOT}/${TEST_SYSTEM}"

REF_LIGAND="${SYSTEM_DIR}/${TEST_SYSTEM}.lig.am1bcc.mol2"
SPHERES="${SYSTEM_DIR}/${TEST_SYSTEM}.rec.clust.close.sph"
GRID_PREFIX="${SYSTEM_DIR}/${TEST_SYSTEM}.rec"
GRID_NRG="${GRID_PREFIX}.nrg"
GRID_BMP="${GRID_PREFIX}.bmp"

# DOCK6.14 executables
DOCK_ROOT="/gpfs/projects/rizzo/iamanor/Collaborations/Carissa/002_Dref_Isolated_6.14_version"
BASE_ROOT="${DOCK_ROOT}/baseline_7d38760d"

BASELINE_DOCK_BIN="${DOCK_ROOT}/gpu_long_builds/dock6.baseline.gpu"
MODULARIZED_DOCK_BIN="${DOCK_ROOT}/gpu_long_builds/dock6.modularized.gpu"

# Input and output
CANONICAL_INPUT="${WORK_ROOT}/dn_generic.in"
OUTPUT_PREFIX="Dref_ouput"

BASELINE_RUN_DIR="${WORK_ROOT}/baseline/system"
MODULARIZED_RUN_DIR="${WORK_ROOT}/modularized/system"

BASELINE_INPUT="${BASELINE_RUN_DIR}/dn_generic.in"
MODULARIZED_INPUT="${MODULARIZED_RUN_DIR}/dn_generic.in"

BASELINE_OUTPUT="${BASELINE_RUN_DIR}/dn_generic.out"
MODULARIZED_OUTPUT="${MODULARIZED_RUN_DIR}/dn_generic.out"

RESULTS_ROOT="${WORK_ROOT}/comparison_results"

BASELINE_LOG="${RESULTS_ROOT}/baseline_run.log"
MODULARIZED_LOG="${RESULTS_ROOT}/modularized_run.log"

COMPARISON_REPORT="${RESULTS_ROOT}/equivalence_report.txt"

# SLURM resources
SLURM_PARTITION="gpu-long"
SLURM_TIME="2-00:00:00"
SLURM_NODES=1
SLURM_NTASKS=1

# Read a DOCK input parameter
parameter_value()
{
    local KEY="$1"
    awk -v key="$KEY" '$1 == key {value=$2; count++} END {if (count == 1) print value}' "$CANONICAL_INPUT"
}

# Validate a required parameter
check_parameter()
{
    local KEY="$1"
    local EXPECTED="$2"
    local ACTUAL

    ACTUAL=$(parameter_value "$KEY")

    if [[ "$ACTUAL" != "$EXPECTED" ]]; then
        echo "ERROR: ${KEY}: expected ${EXPECTED}, found ${ACTUAL}"
        return 1
    fi
}

# Check required files
check_file()
{
    local NAME="$1"
    local FILE="$2"

    if [[ -z "$FILE" || ! -f "$FILE" ]]; then
        echo "ERROR: Missing ${NAME}: ${FILE}"
        return 1
    fi
}

# Validate executables and molecular inputs
validate_configuration()
{
    local ERRORS=0
    local KEY
    local FILE

    [[ -x "$BASELINE_DOCK_BIN" ]] || { echo "ERROR: Baseline executable missing"; ERRORS=$((ERRORS+1)); }
    [[ -x "$MODULARIZED_DOCK_BIN" ]] || { echo "ERROR: Modularized executable missing"; ERRORS=$((ERRORS+1)); }

    check_file "canonical input" "$CANONICAL_INPUT" || return 1
    check_file "1O86 reference ligand" "$REF_LIGAND" || ERRORS=$((ERRORS+1))
    check_file "1O86 receptor spheres" "$SPHERES" || ERRORS=$((ERRORS+1))
    check_file "1O86 grid NRG" "$GRID_NRG" || ERRORS=$((ERRORS+1))
    check_file "1O86 grid BMP" "$GRID_BMP" || ERRORS=$((ERRORS+1))

    check_parameter "conformer_search_type" "denovo" || ERRORS=$((ERRORS+1))
    check_parameter "dn_mode" "sample" || ERRORS=$((ERRORS+1))
    check_parameter "dn_sampling_method" "graph" || ERRORS=$((ERRORS+1))
    check_parameter "dn_user_specified_anchor" "yes" || ERRORS=$((ERRORS+1))
    check_parameter "descriptor_score_primary" "yes" || ERRORS=$((ERRORS+1))
    check_parameter "descriptor_use_hungarian" "yes" || ERRORS=$((ERRORS+1))
    check_parameter "descriptor_use_dynamic_reference" "yes" || ERRORS=$((ERRORS+1))
    check_parameter "dn_write_growth_trees" "yes" || ERRORS=$((ERRORS+1))
    check_parameter "simplex_random_seed" "0" || ERRORS=$((ERRORS+1))
    check_parameter "dn_output_prefix" "$OUTPUT_PREFIX" || ERRORS=$((ERRORS+1))

    check_parameter "receptor_site_file" "$SPHERES" || ERRORS=$((ERRORS+1))
    check_parameter "descriptor_grid_score_grid_prefix" "$GRID_PREFIX" || ERRORS=$((ERRORS+1))
    check_parameter "descriptor_hms_score_ref_filename" "$REF_LIGAND" || ERRORS=$((ERRORS+1))

    for KEY in \
        dn_fraglib_scaffold_file \
        dn_fraglib_linker_file \
        dn_fraglib_sidechain_file \
        dn_fraglib_anchor_file \
        dn_torenv_table \
        vdw_defn_file \
        flex_defn_file \
        flex_drive_file
    do
        FILE=$(parameter_value "$KEY")
        check_file "$KEY" "$FILE" || ERRORS=$((ERRORS+1))
    done

    if (( ERRORS > 0 )); then
        echo "Configuration validation failed with ${ERRORS} error(s)."
        return 1
    fi

    echo "Configuration validation passed for ${TEST_SYSTEM}."
}
