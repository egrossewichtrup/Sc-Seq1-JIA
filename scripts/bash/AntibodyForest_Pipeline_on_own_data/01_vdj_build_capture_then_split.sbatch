#!/bin/bash
#SBATCH --job-name=vdj_build_capture_split
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/03_logs/%x_%j.out
#SBATCH --error=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/03_logs/%x_%j.err

set -euo pipefail

project_root="/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline"
core="/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2025-12_core_data"

gro_multi_root="/hpc/dla_lti/yermanos_group/GRO11679/alignment_CR/cellranger_multi_2025-12-01"

conda_sh="/hpc/dla_lti/yermanos_group/software/venvs4/miniforge3/etc/profile.d/conda.sh"
r_env="/hpc/dla_lti/yermanos_group/software/venvs4/r"

cfg="${project_root}/00_config"
inputs="${project_root}/01_inputs"
results="${project_root}/02_results/vdj_build"
logs="${project_root}/03_logs"
scripts="${project_root}/04_scripts"

mkdir -p "${cfg}" "${inputs}" "${results}" "${logs}" "${scripts}"

export TMPDIR="${TMPDIR:-/scratch/${SLURM_JOB_ID}}"
mkdir -p "${TMPDIR}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

echo "Job started: $(date)"
echo "HOST: $(hostname)"
echo "PROJECT: ${project_root}"
echo "CORE: ${core}"
echo "GRO_MULTI_ROOT: ${gro_multi_root}"
echo "TMPDIR: ${TMPDIR}"
echo "CPUS: ${SLURM_CPUS_PER_TASK}"

captures_tsv="${cfg}/captures_full_vdj_dirs.tsv"
samples_tsv="${cfg}/samples_demux_barcodes.tsv"

# ------------------------------------------------------------
# Create capture mapping if missing
# ------------------------------------------------------------
if [[ ! -s "${captures_tsv}" ]]; then
  echo -e "capture_id\tvdj_full_dir" > "${captures_tsv}"

  for cid in C1 C2 C3 C4 C5 C6 C7 C8 C9 C10; do
    found=""

    if [[ -d "${gro_multi_root}" ]]; then
      found="$(find "${gro_multi_root}" -type d -name vdj_b -path "*/${cid} */*/outs/per_sample_outs/*/vdj_b" 2>/dev/null | head -n 1 || true)"
    fi

    if [[ -z "${found}" ]]; then
      if [[ -d "${core}/${cid}/VDJ" ]]; then
        found="${core}/${cid}/VDJ"
      fi
    fi

    if [[ -z "${found}" ]]; then
      found="__MISSING__"
    fi

    echo -e "${cid}\t${found}" >> "${captures_tsv}"
  done

  echo "Created ${captures_tsv}"
fi

# ------------------------------------------------------------
# Create sample demux mapping if missing
# ------------------------------------------------------------
if [[ ! -s "${samples_tsv}" ]]; then
  cat > "${samples_tsv}" <<EOF
sample_id	patient_id	tissue	side	capture_id	demux_filtered_contig_annotations_csv
HC001_B	HC001	blood	NA	C1	${core}/C1/VDJ/C1_HC001/filtered_contig_annotations.csv
HC002_B	HC002	blood	NA	C1	${core}/C1/VDJ/C1_HC002/filtered_contig_annotations.csv
P001_LK	P001	SF	LK	C2	${core}/C2/VDJ/filtered_contig_annotations.csv
P001_RK	P001	SF	RK	C3	${core}/C3/VDJ/filtered_contig_annotations.csv
P002_LK	P002	SF	LK	C4	${core}/C4/VDJ/filtered_contig_annotations.csv
P002_RK	P002	SF	RK	C5	${core}/C5/VDJ/filtered_contig_annotations.csv
P003_B	P003	blood	NA	C6	${core}/C6/VDJ/C6_P003/filtered_contig_annotations.csv
HC003_B	HC003	blood	NA	C6	${core}/C6/VDJ/C6_HC003/filtered_contig_annotations.csv
P004_LK	P004	SF	LK	C7	${core}/C7/VDJ/filtered_contig_annotations.csv
P005_B	P005	blood	NA	C8	${core}/C8/VDJ/C8_P005/filtered_contig_annotations.csv
P006_B	P006	blood	NA	C8	${core}/C8/VDJ/C8_P006/filtered_contig_annotations.csv
P002_B	P002	blood	NA	C9	${core}/C9/VDJ/C9_P002/filtered_contig_annotations.csv
P007_B	P007	blood	NA	C9	${core}/C9/VDJ/C9_P007/filtered_contig_annotations.csv
P008_B	P008	blood	NA	C9	${core}/C9/VDJ/C9_P008/filtered_contig_annotations.csv
P009_B	P009	blood	NA	C10	${core}/C10/VDJ/C10_P009/filtered_contig_annotations.csv
P010_B	P010	blood	NA	C10	${core}/C10/VDJ/C10_P010/filtered_contig_annotations.csv
EOF
  echo "Created ${samples_tsv}"
fi

# ------------------------------------------------------------
# Fixed missing capture check
# ------------------------------------------------------------
if grep -q $'\t__MISSING__' "${captures_tsv}"; then
  echo "ERROR: Some captures are missing a vdj_b directory in ${captures_tsv}"
  echo "Edit that file and set vdj_full_dir for the missing captures."
  exit 1
fi

# ------------------------------------------------------------
# Activate R and run
# ------------------------------------------------------------
source "${conda_sh}"
conda activate "${r_env}"

r_script="${scripts}/01_vdj_build_capture_then_split.R"

cat > "${r_script}" <<'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Platypus)
})

project_root <- "/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline"
cfg <- file.path(project_root, "00_config")
results_root <- file.path(project_root, "02_results", "vdj_build")

captures_tsv <- file.path(cfg, "captures_full_vdj_dirs.tsv")
samples_tsv  <- file.path(cfg, "samples_demux_barcodes.tsv")

dir.create(results_root, showWarnings = FALSE, recursive = TRUE)

read_tsv_base <- function(path) {
  x <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(x) < 1) stop("Empty TSV: ", path)
  x
}

cap <- read_tsv_base(captures_tsv)
smp <- read_tsv_base(samples_tsv)

must_have <- c(
  "filtered_contig_annotations.csv",
  "filtered_contig.fasta",
  "consensus_annotations.csv",
  "consensus.fasta",
  "clonotypes.csv",
  "concat_ref.fasta"
)

normalize_barcode <- function(x) {
  x <- as.character(x)
  sub("-1$", "", x)
}

assert_full_vdj_dir <- function(d, capture_id) {
  if (!dir.exists(d)) stop("Full VDJ dir not found for ", capture_id, ": ", d)
  missing <- must_have[!file.exists(file.path(d, must_have))]
  if (length(missing) > 0) {
    stop(
      "Full VDJ dir for ", capture_id, " is missing files: ",
      paste(missing, collapse = ", "),
      "\nDir: ", d
    )
  }
}

qc_lines <- list()

for (cid in unique(cap$capture_id)) {
  d <- cap$vdj_full_dir[cap$capture_id == cid][1]
  assert_full_vdj_dir(d, cid)

  cat("\n============================================================\n")
  cat("Running VDJ_build for capture:", cid, "\n")
  cat("Full dir:", d, "\n")

  vdj_cap <- VDJ_build(
    VDJ.sample.list        = c(d),
    remove.divergent.cells = TRUE,
    complete.cells.only    = TRUE,
    trim.germlines         = TRUE
  )

  vdj_cap$capture_id <- cid

  if ("clonotype_id" %in% colnames(vdj_cap)) {
    vdj_cap$clonotype_id_10x_raw <- vdj_cap$clonotype_id
    vdj_cap$clonotype_id <- paste0(cid, "__", vdj_cap$clonotype_id)
  }

  cap_dir <- file.path(results_root, "by_capture", cid)
  dir.create(cap_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(vdj_cap, file.path(cap_dir, paste0("VDJ_build_capture_", cid, ".rds.xz")), compress = "xz")

  if (!("barcode" %in% colnames(vdj_cap))) stop("VDJ_build output has no barcode column for capture ", cid)
  vdj_cap$barcode_key <- normalize_barcode(vdj_cap$barcode)

  smp_sub <- smp[smp$capture_id == cid, , drop = FALSE]
  if (nrow(smp_sub) < 1) next

  for (i in seq_len(nrow(smp_sub))) {
    sid <- smp_sub$sample_id[i]
    pid <- smp_sub$patient_id[i]
    tissue <- smp_sub$tissue[i]
    side <- smp_sub$side[i]
    demux_csv <- smp_sub$demux_filtered_contig_annotations_csv[i]

    if (!file.exists(demux_csv)) stop("Missing demux filtered_contig_annotations.csv for ", sid, ": ", demux_csv)

    demux <- read.csv(demux_csv, stringsAsFactors = FALSE)
    if (!("barcode" %in% colnames(demux))) stop("Demux file has no barcode column: ", demux_csv)

    keys <- unique(normalize_barcode(demux$barcode))
    keep <- vdj_cap$barcode_key %in% keys
    vdj_s <- vdj_cap[keep, , drop = FALSE]

    vdj_s$sample_id <- sid
    vdj_s$patient_id <- pid
    vdj_s$tissue <- tissue
    vdj_s$side <- side

    out_dir <- file.path(results_root, "by_sample", sid)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    out_rds <- file.path(out_dir, paste0("VDJ_build_", sid, ".rds.xz"))
    saveRDS(vdj_s, out_rds, compress = "xz")

    qc_lines[[length(qc_lines) + 1]] <- data.frame(
      capture_id = cid,
      sample_id = sid,
      patient_id = pid,
      tissue = tissue,
      side = side,
      n_rows = nrow(vdj_s),
      n_barcodes_demux = length(keys),
      stringsAsFactors = FALSE
    )

    cat("Saved sample:", sid, "rows:", nrow(vdj_s), "\n")
  }
}

qc <- do.call(rbind, qc_lines)
qc_path <- file.path(results_root, "QC_counts_by_sample.tsv")
write.table(qc, qc_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWrote QC table:", qc_path, "\n")

patients <- unique(smp$patient_id)
for (pid in patients) {
  smp_ids <- smp$sample_id[smp$patient_id == pid]
  paths <- file.path(results_root, "by_sample", smp_ids, paste0("VDJ_build_", smp_ids, ".rds.xz"))
  ok <- file.exists(paths)
  if (any(!ok)) stop("Missing per sample RDS for patient ", pid, ": ", paste(smp_ids[!ok], collapse = ", "))

  lst <- lapply(paths, readRDS)
  vdj_p <- do.call(rbind, lst)

  out_dir <- file.path(results_root, "by_patient", pid)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_rds <- file.path(out_dir, paste0("VDJ_build_patient_", pid, ".rds.xz"))
  saveRDS(vdj_p, out_rds, compress = "xz")

  cat("Saved patient:", pid, "rows:", nrow(vdj_p), "\n")
}

cat("\nDONE\n")
RSCRIPT

chmod +x "${r_script}"

echo "Running R script: ${r_script}"
Rscript "${r_script}"

echo "Job finished: $(date)"
