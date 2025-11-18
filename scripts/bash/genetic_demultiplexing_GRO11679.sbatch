

###### EXAMPLE SCRIPT FOR CAPTURE C2 - REMEBER TOCHANGE DIRETORIES AND K ###########





########## STEP 0 — set project path + make a tmp folder #########
mkdir -p "/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing"
export PROJ="/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing"
mkdir -p "$PROJ/scripts" "$PROJ/demux" "$PROJ/tmp"

########## END STEP 0 #########









######### STEP 1 — create the Slurm job for C2 (K=2) ##########

cat > "$PROJ/scripts/run_souporcell_C2.sbatch" <<'SBATCH'
#!/usr/bin/env bash
#SBATCH --job-name=souporcell_C2
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --chdir=/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL

set -eo pipefail
set -u

PROJ="/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348"

# --- C2 inputs (per-sample outputs) ---
BAM="/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/cellranger_multi_RERUN_BAM/multi_C2/multi_C2_bam/outs/per_sample_outs/multi_C2_bam/count/sample_alignments.bam"
BARCODES_GZ="/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/cellranger_multi_RERUN_BAM/multi_C2/multi_C2_bam/outs/per_sample_outs/multi_C2_bam/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"
REF_FASTA="/hpc/dla_lti/yermanos_group/REF_GEX/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

# --- donors in this capture (change if needed) ---
K=2

# Fresh timestamped outdir
TS="$(date +%Y%m%d_%H%M%S)"
OUTROOT="${PROJ}/demux/genetic_demultiplexing_${TS}"
OUTDIR="${OUTROOT}/C2_souporcell"
mkdir -p "${OUTDIR}"

# Project-local temp (avoid /scratch and small /tmp)
JOBTMP="${PROJ}/tmp/souporcell_${SLURM_JOB_ID}_${TS}"
mkdir -p "$JOBTMP"
export TMPDIR="$JOBTMP"
export TMP="$JOBTMP"
export TEMP="$JOBTMP"
export BCFTOOLS_TMP_DIR="$JOBTMP"

echo "[INFO] Output root: ${OUTROOT}"
echo "[INFO] Souporcell outdir: ${OUTDIR}"
echo "[INFO] Temp dir: ${JOBTMP}"

# Conda activation (disable nounset just here)
set +u
source /hpc/dla_lti/yermanos_group/GRO11679/analysis_20251002_112741/venvs/miniforge3/etc/profile.d/conda.sh
conda activate "${PROJ}/venvs/souporcell"
set -u

# Sanity checks
[ -s "${BAM}" ] || { echo "[ERROR] BAM missing: ${BAM}" >&2; exit 2; }
[ -s "${BARCODES_GZ}" ] || { echo "[ERROR] barcodes.tsv.gz missing: ${BARCODES_GZ}" >&2; exit 2; }
[ -s "${REF_FASTA}" ] || { echo "[ERROR] FASTA missing: ${REF_FASTA}" >&2; exit 2; }

# Show space in project (temp+out live here)
df -h "$PROJ" || true

# Ensure FASTA index
[ -f "${REF_FASTA}.fai" ] || samtools faidx "${REF_FASTA}"

# Prepare barcodes (Souporcell expects uncompressed)
gunzip -c "${BARCODES_GZ}" > "${OUTDIR}/barcodes.tsv"

# Ensure BAM index
[ -f "${BAM}.bai" ] || samtools index "${BAM}"

# Run Souporcell
echo "[INFO] Running Souporcell (threads=${SLURM_CPUS_PER_TASK}, K=${K})"
souporcell_pipeline.py \
  -i "${BAM}" \
  -b "${OUTDIR}/barcodes.tsv" \
  -f "${REF_FASTA}" \
  -t "${SLURM_CPUS_PER_TASK}" \
  -o "${OUTDIR}" \
  -k "${K}"

echo "[DONE] Souporcell C2 finished at $(date)"
echo "[DONE] Outputs in: ${OUTDIR}"
SBATCH

######### END STEP 1 — create the Slurm job for C2 (K=2) ##########






####### STEP 2 — submit the job #######

cd "$PROJ/scripts"
JOBID=$(sbatch run_souporcell_C2.sbatch | awk '{print $4}')
echo "Submitted C2 JobID: $JOBID"

####### END STEP 2 — submit the job #######






####### GENERATE VALIDATION BAR PLOT #######

source /hpc/dla_lti/yermanos_group/GRO11679/analysis_20251002_112741/venvs/miniforge3/etc/profile.d/conda.sh
conda activate /hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/venvs/souporcell


export PROJ="/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing"
mkdir -p "$PROJ/scripts"

cat > "$PROJ/scripts/plot_souporcell_bar_2donors.py" <<'PY'
import argparse, pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--clusters", required=True, help="Path to clusters.tsv from souporcell")
p.add_argument("--out", required=True, help="Output image (e.g. .png)")
p.add_argument("--title", default="Souporcell demultiplexing (K=2)")
p.add_argument("--person1-label", default="person 1")
p.add_argument("--person2-label", default="person 2")
args = p.parse_args()

df = pd.read_csv(args.clusters, sep="\t")

if "assignment" not in df.columns:
    raise SystemExit("clusters.tsv missing 'assignment' column")

assign = df["assignment"].astype(str)
status = df["status"].astype(str) if "status" in df.columns else pd.Series([""]*len(df))

# Doublet if assignment like "0/1" or status contains "doublet"
is_doublet = assign.str.contains("/", regex=False) | status.str.contains("doublet", case=False, na=False)
is_p1 = assign.eq("0")
is_p2 = assign.eq("1")
is_neg = ~(is_doublet | is_p1 | is_p2)

counts = {
    "Doublet": int(is_doublet.sum()),
    args.person1_label: int(is_p1.sum()),
    args.person2_label: int(is_p2.sum()),
    "negative": int(is_neg.sum()),
}
order = ["Doublet", args.person1_label, args.person2_label, "negative"]
heights = [counts[k] for k in order]

plt.figure(figsize=(6,4))
plt.bar(order, heights)
plt.ylabel("Cells")
plt.title(args.title)
plt.tight_layout()
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(args.out, dpi=200)
print("Counts:", counts)
print("Wrote:", args.out)
PY

#Find out directory
LAST_OUT=$(ls -t "$PROJ"/souporcell_C2_*.out | head -1)
echo "OUTLOG: $LAST_OUT"
OUTDIR=$(grep -m1 'Souporcell outdir:' "$LAST_OUT" | awk '{print $4}')
echo "OUTDIR: $OUTDIR"

#Make Barplot
python "$PROJ/scripts/plot_souporcell_bar_2donors.py" \
  --clusters "$OUTDIR/clusters.tsv" \
  --out "$OUTDIR/souporcell_bar_C2.png" \
  --title "Capture C2 demultiplexing (K=2)" \
  --person1-label "person 1" \
  --person2-label "person 2"

ls -lh "$OUTDIR/souporcell_bar_C2.png"

####### END GENERATE VALIDATION BAR PLOT #######