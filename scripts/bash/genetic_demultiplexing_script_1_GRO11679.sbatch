#!/usr/bin/env bash
# Launcher script (run on login node):
#   bash run_souporcell_launcher.sh
#
# This version submits ONE Slurm array job, one task per capture.

set -euo pipefail

#####################################
# === USER SETTINGS (edit here) === #
#####################################

# Top-level demultiplexing root
DEMUX_ROOT="/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing2"

# Project root for genetic demultiplexing with Souporcell
# -> all outputs, logs, tmp & sbatch files will live under here
PROJ="${DEMUX_ROOT}"

# 10x reference FASTA (Souporcell uses this; will index if needed)
REF_FASTA="/hpc/dla_lti/yermanos_group/REF_GEX/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

# Cell Ranger multi root (top-level directory containing the per-capture folders)
CR_RUN_ROOT="/hpc/dla_lti/yermanos_group/GRO11679/alignment_CR/cellranger_multi_2025-12-01"

# Conda (Miniforge) and env for souporcell
MINIFORGE_ROOT="/hpc/dla_lti/yermanos_group/software/venvs4/miniforge3"
CONDA_ENV_PATH="/hpc/dla_lti/yermanos_group/software/venvs4/souporcell"

# >>> Captures + K values (same length, same order) <<<
CAPTURE_IDS=(C1 C6)
K_VALUES=(2  2)

# Optional email notifications (leave empty to disable)
MAIL_USER="${MAIL_USER:-egrossewichtrup@umcutrecht.nl}"

#####################################
# === (Usually no edits below)  === #
#####################################

# Basic sanity: same number of CAPTURE_IDS and K_VALUES
N_JOBS=${#CAPTURE_IDS[@]}
if [[ ${#K_VALUES[@]} -ne ${N_JOBS} ]]; then
  echo "[FATAL] CAPTURE_IDS and K_VALUES must have the same length." >&2
  echo "        CAPTURE_IDS: ${#CAPTURE_IDS[@]}  K_VALUES: ${#K_VALUES[@]}" >&2
  exit 1
fi

ARRAY_MAX=$((N_JOBS - 1))

# Create project subdirs inside demultiplexing root
mkdir -p "${PROJ}" "${PROJ}/demux" "${PROJ}/tmp" "${PROJ}/logs" "${PROJ}/sbatch"

# Check key shared paths once
missing=0
for p in "$REF_FASTA" "$MINIFORGE_ROOT" "$CONDA_ENV_PATH" "$CR_RUN_ROOT"; do
  if [[ ! -e "$p" ]]; then
    echo "[ERROR] Path not found: $p" >&2
    missing=1
  fi
done
if (( missing )); then
  echo "[FATAL] One or more required shared paths are missing – aborting." >&2
  exit 1
fi

TS="$(date +%Y%m%d_%H%M%S)"
SBATCH_FILE="${PROJ}/sbatch/run_souporcell_array_${TS}.sbatch"

echo "[INFO] Writing array sbatch file to: ${SBATCH_FILE}"

cat > "$SBATCH_FILE" <<SBATCH
#!/usr/bin/env bash
#SBATCH --job-name=souporcell_multi
#SBATCH --array=0-${ARRAY_MAX}
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --chdir=${PROJ}
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#MAILTAG

set -euo pipefail

# Values baked in by the launcher
PROJ="${PROJ}"
REF_FASTA="${REF_FASTA}"
CR_RUN_ROOT="${CR_RUN_ROOT}"
MINIFORGE_ROOT="${MINIFORGE_ROOT}"
CONDA_ENV_PATH="${CONDA_ENV_PATH}"

# Per-task capture + K (order must match launcher arrays)
CAPTURE_IDS=(${CAPTURE_IDS[*]})
K_VALUES=(${K_VALUES[*]})

TASK_ID=\${SLURM_ARRAY_TASK_ID}
CAPTURE_ID=\${CAPTURE_IDS[\${TASK_ID}]}
K=\${K_VALUES[\${TASK_ID}]}

echo "[INFO] Starting souporcell for capture \${CAPTURE_ID} (K=\${K}) at: \$(date)"

OUTROOT="\${PROJ}/demux/\${CAPTURE_ID}"
OUTDIR="\${OUTROOT}/\${CAPTURE_ID}_souporcell"
mkdir -p "\${OUTDIR}"

# Job-local temp (one per array task)
JOBTMP="\${PROJ}/tmp/souporcell_\${CAPTURE_ID}_\${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}_${TS}"
mkdir -p "\${JOBTMP}"
export TMPDIR="\${JOBTMP}"
export TMP="\${JOBTMP}"
export TEMP="\${JOBTMP}"
export BCFTOOLS_TMP_DIR="\${JOBTMP}"

echo "[INFO] Output root: \${OUTROOT}"
echo "[INFO] Outdir:      \${OUTDIR}"
echo "[INFO] Temp dir:    \${JOBTMP}"

# Activate conda env
set +u
source "\${MINIFORGE_ROOT}/etc/profile.d/conda.sh"
conda activate "\${CONDA_ENV_PATH}"
set -u

# Tool sanity checks
command -v souporcell_pipeline.py >/dev/null 2>&1 || { echo "[ERROR] souporcell_pipeline.py not found in env"; exit 2; }
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found in env"; exit 2; }

########################################
# Derive BAM + barcodes using 'find'   #
# to match the structure you showed    #
########################################

# 1) Find the capture directory:
#    e.g. "\$CR_RUN_ROOT/C1 ScSeq1-Capture1-HC001-HC002-10-07-2025"
CAPTURE_DIR=\$(find "\${CR_RUN_ROOT}" -maxdepth 1 -type d -name "\${CAPTURE_ID} *" | head -n 1)

if [[ -z "\${CAPTURE_DIR}" ]]; then
  echo "[ERROR] Could not find capture directory for \${CAPTURE_ID} under \${CR_RUN_ROOT}" >&2
  exit 2
fi
echo "[INFO] Using capture dir: \${CAPTURE_DIR}"

# 2) Inside that, find the multi directory:
#    e.g. "\$CAPTURE_DIR/multi_C1_ScSeq1-Capture1-HC001-HC002-10-07-2025"
MULTI_DIR=\$(find "\${CAPTURE_DIR}" -maxdepth 1 -type d -name "multi_\${CAPTURE_ID}_*" | head -n 1)

if [[ -z "\${MULTI_DIR}" ]]; then
  echo "[ERROR] Could not find multi dir for \${CAPTURE_ID} under \${CAPTURE_DIR}" >&2
  exit 2
fi

MULTI_BASENAME=\$(basename "\${MULTI_DIR}")
echo "[INFO] Using multi dir: \${MULTI_DIR} (basename=\${MULTI_BASENAME})"

# 3) Build BAM and barcodes paths like in your example:
#    .../outs/per_sample_outs/<MULTI_BASENAME>/count/sample_alignments.bam
BAM="\${MULTI_DIR}/outs/per_sample_outs/\${MULTI_BASENAME}/count/sample_alignments.bam"
BARCODES_GZ="\${MULTI_DIR}/outs/per_sample_outs/\${MULTI_BASENAME}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"

echo "[INFO] BAM path:      \${BAM}"
echo "[INFO] Barcodes path: \${BARCODES_GZ}"

# Input checks
[[ -s "\${REF_FASTA}" ]] || { echo "[ERROR] FASTA missing: \${REF_FASTA}" >&2; exit 2; }
[[ -s "\${BAM}" ]] || { echo "[ERROR] BAM missing: \${BAM}" >&2; exit 2; }
[[ -s "\${BARCODES_GZ}" ]] || { echo "[ERROR] barcodes.tsv.gz missing: \${BARCODES_GZ}" >&2; exit 2; }

# Disk info for PROJ (temp + output live here)
df -h "\${PROJ}" || true

# Index FASTA and BAM if needed
[[ -f "\${REF_FASTA}.fai" ]] || samtools faidx "\${REF_FASTA}"
[[ -f "\${BAM}.bai" ]] || samtools index "\${BAM}"

# Souporcell expects uncompressed barcodes.tsv
BARCODES_TSV="\${OUTDIR}/barcodes.tsv"
if [[ ! -s "\${BARCODES_TSV}" ]]; then
  gunzip -c "\${BARCODES_GZ}" > "\${BARCODES_TSV}"
fi

echo "[INFO] Running Souporcell: capture=\${CAPTURE_ID}, K=\${K}, threads=\${SLURM_CPUS_PER_TASK}"
souporcell_pipeline.py \\
  -i "\${BAM}" \\
  -b "\${BARCODES_TSV}" \\
  -f "\${REF_FASTA}" \\
  -t "\${SLURM_CPUS_PER_TASK}" \\
  -o "\${OUTDIR}" \\
  -k "\${K}"

echo "[DONE] Souporcell \${CAPTURE_ID} finished at \$(date)"
echo "[DONE] Outputs in: \${OUTDIR}"
SBATCH

# Inject mail-user line if MAIL_USER is set
if [[ -n "$MAIL_USER" ]]; then
  sed -i "s|#MAILTAG|#SBATCH --mail-user=${MAIL_USER}|" "$SBATCH_FILE"
else
  sed -i "s|#MAILTAG||" "$SBATCH_FILE"
fi

chmod +x "$SBATCH_FILE"

echo "[INFO] Submitting array job: $SBATCH_FILE"
JOBID=$(sbatch "$SBATCH_FILE" | awk '{print $4}')
echo "[SUBMITTED] souporcell_multi array — JobID: ${JOBID}"
echo "[TIP] Logs per task will be in ${PROJ}/logs:"
echo "      logs/souporcell_multi_${JOBID}_<task>.out"
echo "      logs/souporcell_multi_${JOBID}_<task>.err"
