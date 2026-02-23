#!/bin/bash
#SBATCH --job-name=platypus_vgm_16
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2025-12_vgm/logs/%x_%j.out
#SBATCH --error=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2025-12_vgm/logs/%x_%j.err

set -euo pipefail

# ----------------------------
# Paths / environment
# ----------------------------
work_dir="/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2025-12_vgm"
inputs_tsv="${work_dir}/resolved_platypus_inputs.tsv"
conda_sh="/hpc/dla_lti/yermanos_group/software/venvs4/miniforge3/etc/profile.d/conda.sh"
r_env="/hpc/dla_lti/yermanos_group/software/venvs4/r"

mkdir -p "${work_dir}/logs"
cd "${work_dir}"

# Use node-local scratch if available
export TMPDIR="${TMPDIR:-/scratch/${SLURM_JOB_ID}}"
mkdir -p "${TMPDIR}"

# Keep BLAS/OMP threads aligned with Slurm allocation (avoid oversubscription)
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

echo "Job started: $(date)"
echo "HOST: $(hostname)"
echo "WORKDIR: ${work_dir}"
echo "TMPDIR: ${TMPDIR}"
echo "CPUS: ${SLURM_CPUS_PER_TASK}"
echo "MEM:  ${SLURM_MEM_PER_NODE:-unknown}"

source "${conda_sh}"
conda activate "${r_env}"

# ----------------------------
# Write R scripts into WORKDIR
# ----------------------------

cat > "${work_dir}/04_make_vgm_16samples.R" <<'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(Seurat)
  library(Platypus)
  library(Matrix)
})

work_dir   <- "/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2025-12_vgm"
inputs_tsv <- file.path(work_dir, "resolved_platypus_inputs.tsv")

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

# Where we will build "outs-like" GEX folders for EVERY sample
gex_outs_root <- file.path(work_dir, "_GEX_outs_like")
dir.create(gex_outs_root, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
barcode_only <- function(x) {
  m <- str_match(x, "([ACGT]+-1)$")[, 2]
  if (any(is.na(m))) {
    stop("Some cell names did not end in a 10x barcode like AAAC...-1.\n",
         "Example bad names:\n", paste(head(x[is.na(m)], 5), collapse = "\n"))
  }
  m
}

get_counts <- function(obj) {
  tryCatch(
    GetAssayData(obj, layer = "counts"),
    error = function(e) GetAssayData(obj, slot = "counts")
  )
}

write_filtered_feature_bc_matrix <- function(counts, outs_dir) {
  mtx_dir <- file.path(outs_dir, "filtered_feature_bc_matrix")
  dir.create(mtx_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Counts matrix must have rownames (genes) and colnames (barcodes).")
  }

  has_gzip <- nzchar(Sys.which("gzip"))

  mtx_path <- file.path(mtx_dir, "matrix.mtx")
  Matrix::writeMM(counts, file = mtx_path)
  if (has_gzip) system2("gzip", c("-f", mtx_path))

  bc_path <- file.path(mtx_dir, "barcodes.tsv")
  writeLines(colnames(counts), con = bc_path)
  if (has_gzip) system2("gzip", c("-f", bc_path))

  feat_path <- file.path(mtx_dir, "features.tsv")
  genes <- rownames(counts)
  feat_df <- data.frame(genes, genes, "Gene Expression")
  write.table(
    feat_df,
    file = feat_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  if (has_gzip) system2("gzip", c("-f", feat_path))

  invisible(outs_dir)
}

assert_10x_mtx_present <- function(mtx_dir) {
  ok <- file.exists(file.path(mtx_dir, "matrix.mtx")) ||
        file.exists(file.path(mtx_dir, "matrix.mtx.gz"))
  ok <- ok && (file.exists(file.path(mtx_dir, "barcodes.tsv")) ||
               file.exists(file.path(mtx_dir, "barcodes.tsv.gz")))
  ok <- ok && (file.exists(file.path(mtx_dir, "features.tsv")) ||
               file.exists(file.path(mtx_dir, "features.tsv.gz")))
  if (!ok) stop("Missing 10x files in: ", mtx_dir)
}

# ----------------------------
# Load resolved inputs
# ----------------------------
tbl <- read_tsv(inputs_tsv, show_col_types = FALSE) %>%
  distinct(sample_id, .keep_all = TRUE)

cat("\nResolved samples:", nrow(tbl), "\n")
if (nrow(tbl) != 16) {
  stop("Expected 16 rows in resolved_platypus_inputs.tsv but got ", nrow(tbl))
}

# ----------------------------
# Build uniform outs-like GEX dirs for ALL samples
# ----------------------------
VDJ.out.directory.list <- vector("list", length = nrow(tbl))
GEX.out.directory.list <- vector("list", length = nrow(tbl))

for (i in seq_len(nrow(tbl))) {
  sid <- tbl$sample_id[i]
  gex_type <- tbl$gex_type[i]

  VDJ.out.directory.list[[i]] <- tbl$vdj_dir[i]

  outs_dir <- file.path(gex_outs_root, sid)
  dir.create(outs_dir, showWarnings = FALSE, recursive = TRUE)

  if (gex_type == "rds") {
    obj <- readRDS(tbl$gex_rds[i])
    if (!("RNA" %in% names(obj@assays))) stop("Sample ", sid, ": no RNA assay found.")
    DefaultAssay(obj) <- "RNA"

    # Make cell names match 10x barcodes
    obj <- RenameCells(obj, new.names = barcode_only(Cells(obj)))

    counts <- as(get_counts(obj), "dgCMatrix")
    write_filtered_feature_bc_matrix(counts, outs_dir)
    assert_10x_mtx_present(file.path(outs_dir, "filtered_feature_bc_matrix"))

  } else if (gex_type == "mtx") {
    # Create a "filtered_feature_bc_matrix" link pointing to the real matrix folder
    real_mtx <- tbl$gex_mtx[i]
    link_path <- file.path(outs_dir, "filtered_feature_bc_matrix")

    # Remove old link/dir if present (safe)
    if (file.exists(link_path) || dir.exists(link_path)) {
      unlink(link_path, recursive = TRUE, force = TRUE)
    }

    ok <- file.symlink(from = real_mtx, to = link_path)
    if (!isTRUE(ok)) stop("Failed to symlink ", link_path, " -> ", real_mtx)

    assert_10x_mtx_present(link_path)

  } else {
    stop("Sample ", sid, ": unexpected gex_type=", gex_type)
  }

  GEX.out.directory.list[[i]] <- outs_dir
}

cat("\n--- PAIRING ORDER (index, VDJ, GEX outs-like) ---\n")
for (i in seq_along(VDJ.out.directory.list)) {
  cat(sprintf("%02d VDJ=%s\n   GEX=%s\n",
              i, VDJ.out.directory.list[[i]], GEX.out.directory.list[[i]]))
}

# ----------------------------
# Vignette-style run (full 16)
# ----------------------------
vgm <- VDJ_GEX_matrix(
  VDJ.out.directory.list = VDJ.out.directory.list,
  GEX.out.directory.list = GEX.out.directory.list,
  GEX.integrate = TRUE,
  VDJ.combine = TRUE,
  integrate.GEX.to.VDJ = TRUE,
  integrate.VDJ.to.GEX = TRUE,
  exclude.GEX.not.in.VDJ = FALSE,
  filter.overlapping.barcodes.GEX = TRUE,
  filter.overlapping.barcodes.VDJ = TRUE,
  get.VDJ.stats = TRUE,
  parallel.processing = "none",
  trim.and.align = FALSE,
  group.id = seq_along(VDJ.out.directory.list)
)

out_rds <- file.path(work_dir, "VGM_16samples_vignette_style.rds")
saveRDS(vgm, out_rds)

cat("\nSUCCESS. Saved:\n", out_rds, "\n")
cat("\nvgm[[1]] rows:", nrow(vgm[[1]]), "\n")
cat("vgm[[2]] cells:", ncol(vgm[[2]]), "\n")
RSCRIPT

chmod +x "${work_dir}/04_make_vgm_16samples.R"

cat > "${work_dir}/05_extract_vdj_only_from_vgm.R" <<'RSCRIPT'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 05_extract_vdj_only_from_vgm.R <VGM_RDS> <OUT_VDJ_RDS> <OUT_MINI_VGM_RDS>")
}
vgm_rds <- args[[1]]
out_vdj <- args[[2]]
out_mini <- args[[3]]

cat("Reading VGM:", vgm_rds, "\n")
v <- readRDS(vgm_rds)

if (!is.list(v) || length(v) < 2) stop("Input does not look like a Platypus VGM list.")

vdj <- v[[1]]
cat("VDJ object class:", paste(class(vdj), collapse=","), "\n")
cat("VDJ rows:", if (!is.null(nrow(vdj))) nrow(vdj) else NA, "\n")
cat("VDJ cols:", if (!is.null(ncol(vdj))) ncol(vdj) else NA, "\n")

# Save JUST the VDJ master dataframe
saveRDS(vdj, out_vdj, compress = "xz")
cat("Saved VDJ-only RDS:", out_vdj, "\n")

# Save a "mini VGM" without the large GEX Seurat object
mini <- list(
  VDJ = v[[1]],
  stats = if (length(v) >= 3) v[[3]] else NULL,
  params = if (length(v) >= 4) v[[4]] else NULL,
  session = if (length(v) >= 5) v[[5]] else NULL
)
saveRDS(mini, out_mini, compress = "xz")
cat("Saved mini-VGM (no GEX):", out_mini, "\n")

cat("Done.\n")
RSCRIPT

chmod +x "${work_dir}/05_extract_vdj_only_from_vgm.R"

# ----------------------------
# Run step 1: create full VGM
# ----------------------------
echo "=== Step 1/2: Creating full 16-sample VGM ==="
Rscript "${work_dir}/04_make_vgm_16samples.R"

vgm_rds="${work_dir}/VGM_16samples_vignette_style.rds"
if [[ ! -s "${vgm_rds}" ]]; then
  echo "ERROR: Expected VGM not found or empty: ${vgm_rds}"
  exit 1
fi
echo "OK: Created ${vgm_rds}"
ls -lh "${vgm_rds}"

# ----------------------------
# Run step 2: extract VDJ-only + mini VGM
# ----------------------------
echo "=== Step 2/2: Extracting VDJ-only and mini-VGM ==="
out_vdj="${work_dir}/VDJ_16samples_only.rds.xz"
out_mini="${work_dir}/VGM_16samples_mini_noGEX.rds.xz"

Rscript "${work_dir}/05_extract_vdj_only_from_vgm.R" "${vgm_rds}" "${out_vdj}" "${out_mini}"

echo "Outputs:"
ls -lh "${out_vdj}" "${out_mini}"

echo "Job finished: $(date)"

