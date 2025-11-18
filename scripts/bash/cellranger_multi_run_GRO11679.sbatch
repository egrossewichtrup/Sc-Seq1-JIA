#!/usr/bin/env bash
#SBATCH --job-name=cr_multi_array
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-5

set -euo pipefail

# --- 1) PROJECT ROOT + DATED RUN FOLDER ---
# Your fixed analysis directory:
PROJ="/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348"
RUN_DATE="$(date +%F)"              # e.g., 2025-10-30
RUN_ROOT="${PROJ}/cellranger_multi_${RUN_DATE}"   # e.g., .../cellranger_multi_2025-10-30

mkdir -p "${RUN_ROOT}" "${PROJ}/logs"

# --- 2) BINARIES / REFERENCES / FASTQ LOCATIONS ---
CR_BIN="/hpc/dla_lti/yermanos_group/software/cellranger-9.0.1/bin"
REF_GEX="/hpc/dla_lti/yermanos_group/REF_GEX/refdata-gex-GRCh38-2024-A"
REF_VDJ="/hpc/dla_lti/yermanos_group/REF_VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
FEATURE_REF="${PROJ}/feature_ref/featurebarcode.csv"

FASTQ_A="/hpc/dla_lti/yermanos_group/GRO11679/fastq/A235V25LT3"
FASTQ_B="/hpc/dla_lti/yermanos_group/GRO11679/fastq/B22WFC2LT3"

# --- 3) CAPTURE LIST (code + human-readable stub) ---
# Outputs will be named with the FULL label (e.g., "C1 ScSeq1-Capture1-HC001-HC002-10-07-2025")
CAPS=(
  "C1 ScSeq1-Capture1-HC001-HC002-10-07-2025"
  "C6 ScSeq1-Capture6-P003-HC003-21-07-2025"
  "C8 ScSeq1-Capture8-P005-P006-24-07-2025"
  "C9 ScSeq1-Capture9-P002blood-P007-P008-30-07-2025"
  "C10 ScSeq1-Capture10-P009-P010-31-07-2025"
)

# --- 4) PICK THE CURRENT CAPTURE FROM THE ARRAY INDEX ---
row="${CAPS[$((SLURM_ARRAY_TASK_ID-1))]}"

# Split into: code (first token) + stub (everything after first space)
code="${row%% *}"
stub="${row#* }"

full_label="${code} ${stub}"                    # e.g., "C1 ScSeq1-Capture1-..."
safe_label="${full_label// /_}"                 # path-safe: "C1_ScSeq1-Capture1-..."

# --- 5) WRITE A PER-CAPTURE multi_config CSV ---
cfg="${RUN_ROOT}/multi_config_${safe_label}.csv"
cat > "${cfg}" <<CFG
[gene-expression]
reference,${REF_GEX}
create-bam,true

[vdj]
reference,${REF_VDJ}

[feature]
reference,${FEATURE_REF}

[libraries]
fastq_id,fastqs,feature_types
${stub}_A235V25LT3,${FASTQ_A},Gene Expression
${stub}_B22WFC2LT3,${FASTQ_B},Gene Expression
${stub}-CSP_A235V25LT3,${FASTQ_A},Antibody Capture
${stub}-CSP_B22WFC2LT3,${FASTQ_B},Antibody Capture
${stub}-BCR_A235V25LT3,${FASTQ_A},VDJ-B
${stub}-BCR_B22WFC2LT3,${FASTQ_B},VDJ-B
CFG

# --- 6) OUT DIR + RUN ---
outdir="${RUN_ROOT}/${full_label}"   # NOTE: this keeps spaces in the folder name (quoted)
mkdir -p "${outdir}"
cd "${outdir}"

"${CR_BIN}/cellranger" multi \
  --id="multi_${safe_label}" \
  --csv="${cfg}"


# Submit with sbatch /hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/scripts/cellranger_multi_run_GRO11679.sbatch
