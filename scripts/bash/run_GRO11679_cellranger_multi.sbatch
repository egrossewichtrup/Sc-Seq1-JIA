#!/usr/bin/env bash
#SBATCH --job-name=GRO11679_cr_multi
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=1-5%2
#SBATCH --chdir=/hpc/dla_lti/yermanos_group/GRO11679/2026-06_alignment_CR
#SBATCH --output=/hpc/dla_lti/yermanos_group/GRO11679/2026-06_alignment_CR/slurm_%x_%A_%a.out
#SBATCH --error=/hpc/dla_lti/yermanos_group/GRO11679/2026-06_alignment_CR/slurm_%x_%A_%a.err

set -euo pipefail

###############################################################################
# 1. PROJECT LOCATIONS
###############################################################################

BASE="/hpc/dla_lti/yermanos_group/GRO11679/2026-06_alignment_CR"

CONFIG_DIR="${BASE}/configs"
OUT_PARENT="${BASE}/cellranger_multi_outputs"
LOG_DIR="${BASE}/logs"
TMP_ROOT="${BASE}/tmp"

mkdir -p \
    "${CONFIG_DIR}" \
    "${OUT_PARENT}" \
    "${LOG_DIR}" \
    "${TMP_ROOT}"

###############################################################################
# 2. CELL RANGER, REFERENCES, AND FASTQ DIRECTORIES
###############################################################################

CR_BIN="/hpc/dla_lti/yermanos_group/software/cellranger-9.0.1/bin"
CELLRANGER="${CR_BIN}/cellranger"

REF_GEX="/hpc/dla_lti/yermanos_group/REF_GEX/refdata-gex-GRCh38-2024-A"
REF_VDJ="/hpc/dla_lti/yermanos_group/REF_VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"

FASTQ_A="/hpc/dla_lti/yermanos_group/GRO11679/fastq/A235V25LT3"
FASTQ_B="/hpc/dla_lti/yermanos_group/GRO11679/fastq/B22WFC2LT3"

# Change this only when intentionally starting a new complete rerun.
RUN_TAG="20260615"

###############################################################################
# 3. CAPTURE DEFINITIONS
###############################################################################

CODES=(
    "C1"
    "C6"
    "C8"
    "C9"
    "C10"
)

STUBS=(
    "ScSeq1-Capture1-HC001-HC002-10-07-2025"
    "ScSeq1-Capture6-P003-HC003-21-07-2025"
    "ScSeq1-Capture8-P005-P006-24-07-2025"
    "ScSeq1-Capture9-P002blood-P007-P008-30-07-2025"
    "ScSeq1-Capture10-P009-P010-31-07-2025"
)

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: Submit this script with sbatch as a Slurm array job."
    exit 1
fi

INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if (( INDEX < 0 || INDEX >= ${#CODES[@]} )); then
    echo "ERROR: Invalid SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}."
    exit 1
fi

CODE="${CODES[$INDEX]}"
STUB="${STUBS[$INDEX]}"

RUN_ID="GRO11679_${CODE}_multi_${RUN_TAG}"
RUN_DIR="${OUT_PARENT}/${RUN_ID}"

###############################################################################
# 4. SAMPLE TO HASHTAG MAPPING
###############################################################################

case "${CODE}" in
    C1)
        SAMPLE_ROWS=$'HC001,HTO1\nHC002,HTO2'
        ;;

    C6)
        SAMPLE_ROWS=$'HC003,HTO1\nP003,HTO2'
        ;;

    C8)
        SAMPLE_ROWS=$'P005,HTO3\nP006,HTO4'
        ;;

    C9)
        SAMPLE_ROWS=$'P002,HTO1\nP007,HTO2\nP008,HTO3'
        ;;

    C10)
        SAMPLE_ROWS=$'P009,HTO3\nP010,HTO4'
        ;;

    *)
        echo "ERROR: No hashtag mapping defined for ${CODE}."
        exit 1
        ;;
esac

###############################################################################
# 5. GENERAL INPUT VALIDATION
###############################################################################

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ -x "${CELLRANGER}" ]] || \
    fail "Cell Ranger is not executable: ${CELLRANGER}"

[[ -d "${REF_GEX}" ]] || \
    fail "GEX reference not found: ${REF_GEX}"

[[ -d "${REF_VDJ}" ]] || \
    fail "VDJ reference not found: ${REF_VDJ}"

[[ -d "${FASTQ_A}" ]] || \
    fail "FASTQ directory not found: ${FASTQ_A}"

[[ -d "${FASTQ_B}" ]] || \
    fail "FASTQ directory not found: ${FASTQ_B}"

###############################################################################
# 6. DEFINE AND VALIDATE FASTQ PREFIXES
###############################################################################

GEX_A="${STUB}_A235V25LT3"
GEX_B="${STUB}_B22WFC2LT3"

BCR_A="${STUB}-BCR_A235V25LT3"
BCR_B="${STUB}-BCR_B22WFC2LT3"

CSP_A="${STUB}-CSP_A235V25LT3"
CSP_B="${STUB}-CSP_B22WFC2LT3"

check_fastq_prefix() {
    local fastq_dir="$1"
    local prefix="$2"
    local r1_file
    local r2_file

    r1_file="$(
        find "${fastq_dir}" \
            -type f \
            -name "${prefix}_S*_R1_001.fastq.gz" \
            -print \
            -quit
    )"

    r2_file="$(
        find "${fastq_dir}" \
            -type f \
            -name "${prefix}_S*_R2_001.fastq.gz" \
            -print \
            -quit
    )"

    [[ -n "${r1_file}" ]] || \
        fail "No R1 FASTQ found for ${prefix} in ${fastq_dir}"

    [[ -n "${r2_file}" ]] || \
        fail "No R2 FASTQ found for ${prefix} in ${fastq_dir}"

    echo "Found FASTQs for: ${prefix}"
}

check_fastq_prefix "${FASTQ_A}" "${GEX_A}"
check_fastq_prefix "${FASTQ_B}" "${GEX_B}"

check_fastq_prefix "${FASTQ_A}" "${BCR_A}"
check_fastq_prefix "${FASTQ_B}" "${BCR_B}"

check_fastq_prefix "${FASTQ_A}" "${CSP_A}"
check_fastq_prefix "${FASTQ_B}" "${CSP_B}"

###############################################################################
# 7. PROTECT EXISTING OUTPUTS
###############################################################################

if [[ -e "${RUN_DIR}" ]]; then
    fail "Output already exists: ${RUN_DIR}. Change RUN_TAG or move the existing directory."
fi

###############################################################################
# 8. WRITE THE FEATURE REFERENCE
###############################################################################

FEATURE_REF="${CONFIG_DIR}/feature_reference_${CODE}_${RUN_TAG}.csv"

cat > "${FEATURE_REF}" <<'FEATURE_EOF'
id,name,read,pattern,sequence,feature_type
HTO1,HTO1,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
HTO2,HTO2,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
HTO3,HTO3,R2,5PNNNNNNNNNN(BC),TTCCGCCTCTCTTTG,Antibody Capture
HTO4,HTO4,R2,5PNNNNNNNNNN(BC),AGTAAGTTCAGCGTA,Antibody Capture
FEATURE_EOF

###############################################################################
# 9. WRITE THE CAPTURE SPECIFIC MULTI CONFIG CSV
###############################################################################

MULTI_CSV="${CONFIG_DIR}/multi_config_${CODE}_${RUN_TAG}.csv"

cat > "${MULTI_CSV}" <<CONFIG_EOF
[gene-expression]
reference,${REF_GEX}
create-bam,true

[vdj]
reference,${REF_VDJ}

[feature]
reference,${FEATURE_REF}

[libraries]
fastq_id,fastqs,feature_types
${GEX_A},${FASTQ_A},Gene Expression
${GEX_B},${FASTQ_B},Gene Expression
${BCR_A},${FASTQ_A},VDJ-B
${BCR_B},${FASTQ_B},VDJ-B
${CSP_A},${FASTQ_A},Antibody Capture
${CSP_B},${FASTQ_B},Antibody Capture

[samples]
sample_id,hashtag_ids
${SAMPLE_ROWS}
CONFIG_EOF

###############################################################################
# 10. CELL RANGER RESOURCE SETTINGS
###############################################################################

LOCALCORES="${SLURM_CPUS_PER_TASK:-16}"

# Slurm gives the job 128 GB. Cell Ranger receives 110 GB so that memory
# remains available for the shell, Slurm processes, tee, and filesystem work.
LOCALMEM_GB=110

ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

export TMPDIR="${TMP_ROOT}/${RUN_ID}_${ARRAY_JOB_ID}_${TASK_ID}"
mkdir -p "${TMPDIR}"

###############################################################################
# 11. REMOVE TEMPORARY FILES ONLY AFTER SUCCESS
###############################################################################

SUCCESS=0

cleanup() {
    if [[ "${SUCCESS}" -eq 1 ]]; then
        echo "Removing successful task temporary directory:"
        echo "${TMPDIR}"
        rm -rf "${TMPDIR}"
    else
        echo "Temporary directory retained for troubleshooting:"
        echo "${TMPDIR}"
    fi
}

trap cleanup EXIT

###############################################################################
# 12. PRINT RUN INFORMATION
###############################################################################

echo "=================================================================="
echo "GRO11679 CELL RANGER MULTI"
echo "=================================================================="
echo "Date:              $(date)"
echo "Host:              $(hostname)"
echo "Slurm job:         ${SLURM_JOB_ID:-not_available}"
echo "Array task:        ${SLURM_ARRAY_TASK_ID}"
echo "Capture:           ${CODE}"
echo "FASTQ stub:        ${STUB}"
echo "Cell Ranger:       ${CELLRANGER}"
echo "Run ID:            ${RUN_ID}"
echo "Output parent:     ${OUT_PARENT}"
echo "Config CSV:        ${MULTI_CSV}"
echo "Feature reference: ${FEATURE_REF}"
echo "TMPDIR:            ${TMPDIR}"
echo "Local cores:       ${LOCALCORES}"
echo "Local memory GB:   ${LOCALMEM_GB}"
echo

"${CELLRANGER}" --version

echo
echo "Disk space:"
df -h \
    "${FASTQ_A}" \
    "${FASTQ_B}" \
    "${OUT_PARENT}" \
    "${TMPDIR}" || true

echo
echo "Feature reference:"
cat "${FEATURE_REF}"

echo
echo "Multi config:"
cat "${MULTI_CSV}"

###############################################################################
# 13. RUN CELL RANGER MULTI
###############################################################################

CONSOLE_LOG="${LOG_DIR}/${RUN_ID}.console.log"

cd "${OUT_PARENT}"

echo
echo "Starting Cell Ranger multi for ${CODE}..."

set +e

"${CELLRANGER}" multi \
    --id="${RUN_ID}" \
    --csv="${MULTI_CSV}" \
    --localcores="${LOCALCORES}" \
    --localmem="${LOCALMEM_GB}" \
    --disable-ui \
    2>&1 | tee "${CONSOLE_LOG}"

PIPE_CODES=("${PIPESTATUS[@]}")

set -e

CR_EXIT="${PIPE_CODES[0]}"
TEE_EXIT="${PIPE_CODES[1]}"

echo
echo "Cell Ranger exit code: ${CR_EXIT}"
echo "tee exit code:         ${TEE_EXIT}"

###############################################################################
# 14. FAILURE REPORTING
###############################################################################

if [[ "${CR_EXIT}" -ne 0 || "${TEE_EXIT}" -ne 0 ]]; then
    echo
    echo "Cell Ranger or console log writing failed."
    echo
    echo "Tail of the main Cell Ranger log:"

    tail -n 200 "${RUN_DIR}/_log" 2>/dev/null || true

    echo
    echo "Potential stage error files:"

    find "${RUN_DIR}" \
        \( -name _errors -o -name errors -o -name stderr \) \
        -print \
        2>/dev/null || true

    echo
    echo "Full console log:"
    echo "${CONSOLE_LOG}"

    exit 1
fi

###############################################################################
# 15. SUCCESS SUMMARY
###############################################################################

SUCCESS=1

echo
echo "SUCCESS"
echo "Capture: ${CODE}"
echo
echo "Main output directory:"
echo "${RUN_DIR}/outs"
echo
echo "Expected sample output directory:"
echo "${RUN_DIR}/outs/per_sample_outs"
echo
echo "Top level output structure:"

find "${RUN_DIR}/outs" \
    -maxdepth 3 \
    -print \
    | sort \
    | head -250 || true
