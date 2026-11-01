#!/bin/bash

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_ROOT="${CONFIG_DIR}"

SYSTEM_LIST="${WORK_ROOT}/LISTOFSYSTEMS.txt"
SYSTEM_ROOT="/gpfs/projects/rizzo/iamanor/Systems_and_Library_Files/001_Systems_Files/002_DRef_Systems"
ANCHOR_ROOT="/gpfs/projects/rizzo/iamanor/DOCK6_Development/Dynamics_Referencing/for_isaac/Dynamic_Reference_for_Isaac/anchors"
N_anc=5

CARISSA_ROOT="/gpfs/projects/rizzo/iamanor/Collaborations/Carissa"
INHERITED_SRC="${CARISSA_ROOT}/001_DRef_6.14_version/dock6_dev"
ISOLATED_SRC="${CARISSA_ROOT}/002_Dref_Isolated_6.14_version/dock6_dev"
INHERITED_DOCK_BIN="${INHERITED_SRC}/bin/dref_inherited"
ISOLATED_DOCK_BIN="${ISOLATED_SRC}/bin/dref_isolated"

# The existing validated input supplies all non-system DOCK parameters.
TEMPLATE_INPUT="${CARISSA_ROOT}/003_Equivalence_Testing/DRef-Equivalence-Testing/dn_generic.in"
RUN_ROOT="${WORK_ROOT}/runs"
RESULTS_ROOT="${WORK_ROOT}/comparison_results"

SLURM_PARTITION="rn-long-40core"
SLURM_TIME="2-00:00:00"
SLURM_CPUS=40
MAX_CONCURRENT=40
