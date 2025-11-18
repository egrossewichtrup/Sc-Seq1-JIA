#!/usr/bin/env bash
#SBATCH --job-name=multi_11735
#SBATCH --partition=cpu
#SBATCH --account=dla_lti
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/hpc/dla_lti/yermanos_group/GRO11735/slurm/%x_%A_%a.out
#SBATCH --error=/hpc/dla_lti/yermanos_group/GRO11735/slurm/%x_%A_%a.err

set -euo pipefail
set -x

CR="/hpc/dla_lti/yermanos_group/software/cellranger-9.0.1/bin/cellranger"
ROOT="/hpc/dla_lti/yermanos_group/GRO11735"
FASTQ_ROOT="$ROOT/fastq"
PAIRS="$FASTQ_ROOT/pairs.tsv"

# Date+time stamped output root
RUN_TS="$(date +%F_%H%M)"
OUT_BASE="$ROOT/cellranger_${RUN_TS}"
mkdir -p "$OUT_BASE/logs"

# This array script expects you to submit with --array=1-N (N = lines in pairs.tsv)
LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PAIRS" || true)"
echo "LINE='$LINE'"
if [ -z "$LINE" ]; then
  echo "ERROR: No line found in $PAIRS for index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

GEX_PREFIX="$(printf '%s\n' "$LINE" | awk '{print $1}')"
VDJ_PREFIX="$(printf '%s\n' "$LINE" | awk '{print $2}')"
echo "Selected GEX_PREFIX=$GEX_PREFIX"
echo "Selected VDJ_PREFIX=$VDJ_PREFIX"

# Verify FASTQs exist for lanes 5–8
for L in 5 6 7 8; do
  ls "$FASTQ_ROOT/${GEX_PREFIX}"_*_L00${L}_*.fastq.gz | head -n 1 >/dev/null
  ls "$FASTQ_ROOT/${VDJ_PREFIX}"_*_L00${L}_*.fastq.gz | head -n 1 >/dev/null
done

# Unique per-sample ID (≤64 chars)
HASH="$(printf '%s' "${GEX_PREFIX}_L5-8" | md5sum | awk '{print $1}' | cut -c1-10)"
ID="GRO11735_${HASH}_L5-8"
DESC="GEX=${GEX_PREFIX};VDJ=${VDJ_PREFIX};lanes=5,6,7,8"

# Clean any stale same-ID artifacts (inside this OUT_BASE only)
rm -rf "$OUT_BASE/$ID" "$OUT_BASE/__${ID}.mro"

# Build CSV inline (lanes 5,6,7,8; VDJ-B)
CSV="$OUT_BASE/${ID}.multi.csv"
cat > "$CSV" <<CSV
[gene-expression]
reference,/hpc/dla_lti/yermanos_group/REF_GEX/refdata-gex-GRCh38-2024-A
create-bam,true

[vdj]
reference,/hpc/dla_lti/yermanos_group/REF_VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

[libraries]
fastq_id,fastqs,feature_types,lanes
${GEX_PREFIX},$FASTQ_ROOT,Gene Expression,5|6|7|8
${VDJ_PREFIX},$FASTQ_ROOT,VDJ-B,5|6|7|8
CSV

# Sanity clean (no 'chain' column; enforce VDJ-B)
sed -i '/^[[:space:]]*chain[[:space:]]*,/d' "$CSV"

echo "=== CSV FINAL ==="; sed -n '1,200p' "$CSV"

# Run
export TMPDIR="/scratch/${SLURM_JOB_ID}"
cd "$OUT_BASE"
echo OUT_BASE="$OUT_BASE"
echo ID="$ID"
echo CSV="$CSV"
"$CR" --version
"$CR" multi \
  --id="$ID" \
  --description="$DESC" \
  --csv="$CSV" \
  --jobmode=local \
  --localcores="$SLURM_CPUS_PER_TASK" \
  --localmem=90 \
  --disable-ui

# Friendly summary
echo "=== Outputs (key) ==="
echo "- Web summary: $OUT_BASE/$ID/outs/per_sample_outs/$ID/web_summary.html"
echo "- GEX H5:      $OUT_BASE/$ID/outs/per_sample_outs/$ID/count/sample_filtered_feature_bc_matrix.h5"
echo "- VDJ-B CSV:   $OUT_BASE/$ID/outs/per_sample_outs/$ID/vdj_b/filtered_contig_annotations.csv"
