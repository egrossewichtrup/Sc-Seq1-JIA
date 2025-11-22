#!/usr/bin/env bash
# Launcher script (run on login node):
#   bash run_souporcell_launcher.sh
#
# You can override CAPTURE_ID and K at call time, e.g.:
#   CAPTURE_ID=C2 K=3 bash run_souporcell_launcher.sh
# Run time is multiple hours

set -euo pipefail

#####################################
# === USER SETTINGS (edit here) === #
#####################################

# Top-level demultiplexing root
DEMUX_ROOT="/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing"

# Project root for genetic demultiplexing with Souporcell
# -> all outputs, logs, tmp & sbatch files will live under here
PROJ="${DEMUX_ROOT}"

# 10x reference FASTA (Souporcell uses this; will index if needed)
REF_FASTA="/hpc/dla_lti/yermanos_group/REF_GEX/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

# Cell Ranger multi root (with per-capture multi_<CAPTURE_ID>_bam dirs)
CR_RUN_ROOT="/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/cellranger_multi_RERUN_BAM"

# Conda (Miniforge) and env for souporcell
MINIFORGE_ROOT="/hpc/dla_lti/yermanos_group/software/venvs4/miniforge3"
CONDA_ENV_PATH="/hpc/dla_lti/yermanos_group/software/venvs4/souporcell"

# Capture ID + number of donors
CAPTURE_ID="${CAPTURE_ID:-C1}"
K="${K:-2}"

# Optional email notifications (leave empty to disable)
MAIL_USER="${MAIL_USER:-egrossewichtrup@umcutrecht.nl}"

#####################################
# === (Usually no edits below)  === #
#####################################

# Create project subdirs **inside demultiplexing root**
mkdir -p "${PROJ}" "${PROJ}/demux" "${PROJ}/tmp" "${PROJ}/logs" "${PROJ}/sbatch"

# If BAM/BARCODES_GZ not provided, derive them from CR_RUN_ROOT + CAPTURE_ID
if [[ -z "${BAM:-}" ]]; then
  BAM="${CR_RUN_ROOT}/multi_${CAPTURE_ID}/multi_${CAPTURE_ID}_bam/outs/per_sample_outs/multi_${CAPTURE_ID}_bam/count/sample_alignments.bam"
fi
if [[ -z "${BARCODES_GZ:-}" ]]; then
  BARCODES_GZ="${CR_RUN_ROOT}/multi_${CAPTURE_ID}/multi_${CAPTURE_ID}_bam/outs/per_sample_outs/multi_${CAPTURE_ID}_bam/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"
fi

# Basic checks (fail early if something is missing)
missing=0
for p in "$REF_FASTA" "$BAM" "$BARCODES_GZ" "$MINIFORGE_ROOT"; do
  if [[ ! -e "$p" ]]; then
    echo "[ERROR] Path not found: $p" >&2
    missing=1
  fi
done
if (( missing )); then
  echo "[FATAL] One or more required paths are missing – aborting." >&2
  exit 1
fi

# Timestamp and sbatch filename under PROJ/sbatch
TS="$(date +%Y%m%d_%H%M%S)"
SBATCH_FILE="${PROJ}/sbatch/run_souporcell_${CAPTURE_ID}_${TS}.sbatch"

echo "[INFO] Writing sbatch file to: ${SBATCH_FILE}"

cat > "$SBATCH_FILE" <<SBATCH
#!/usr/bin/env bash
#SBATCH --job-name=souporcell_${CAPTURE_ID}
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --chdir=${PROJ}
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#MAILTAG

set -euo pipefail

echo "[INFO] Starting souporcell for capture ${CAPTURE_ID} at: \$(date)"

# Values baked in by the launcher
PROJ="${PROJ}"
REF_FASTA="${REF_FASTA}"
BAM="${BAM}"
BARCODES_GZ="${BARCODES_GZ}"
CAPTURE_ID="${CAPTURE_ID}"
K="${K}"
MINIFORGE_ROOT="${MINIFORGE_ROOT}"
CONDA_ENV_PATH="${CONDA_ENV_PATH}"

# Timestamp from launcher
OUTROOT="\${PROJ}/demux/genetic_demultiplexing_${TS}"
OUTDIR="\${OUTROOT}/\${CAPTURE_ID}_souporcell"
mkdir -p "\${OUTDIR}"

# Job-local temp (uses child job's SLURM_JOB_ID) – also inside PROJ
JOBTMP="\${PROJ}/tmp/souporcell_\${CAPTURE_ID}_\${SLURM_JOB_ID}_${TS}"
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

echo "[INFO] Running Souporcell: K=\${K}, threads=\${SLURM_CPUS_PER_TASK}"
souporcell_pipeline.py \
  -i "\${BAM}" \
  -b "\${BARCODES_TSV}" \
  -f "\${REF_FASTA}" \
  -t "\${SLURM_CPUS_PER_TASK}" \
  -o "\${OUTDIR}" \
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

echo "[INFO] Submitting sbatch file: $SBATCH_FILE"
JOBID=$(sbatch "$SBATCH_FILE" | awk '{print $4}')
echo "[SUBMITTED] ${CAPTURE_ID} — JobID: \${JOBID}"
echo "[TIP] Child logs will be in ${PROJ}/logs:"
echo "      logs/souporcell_${CAPTURE_ID}_\${JOBID}.out"
echo "      logs/souporcell_${CAPTURE_ID}_\${JOBID}.err"

