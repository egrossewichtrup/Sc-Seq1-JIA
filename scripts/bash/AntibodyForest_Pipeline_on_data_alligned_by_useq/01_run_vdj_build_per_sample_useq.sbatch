#!/bin/bash
#SBATCH --job-name=vdj_build_per_sample_useq_trimmed
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-04_AF_Pipeline_with_useq_aligned_data/03_logs/%x_%j.out
#SBATCH --error=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-04_AF_Pipeline_with_useq_aligned_data/03_logs/%x_%j.err

set -euo pipefail

project_root="/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-04_AF_Pipeline_with_useq_aligned_data"

gro11679_root="/hpc/dla_lti/yermanos_group/GRO11679/aligned_and_demultiplexed_by_useq_GRO11679/unpacked_useq/GRO11679_GRCh38_cellranger"
gro11735_root="/hpc/dla_lti/yermanos_group/GRO11735/aligned_and_demultiplexed_by_useq_GRO11735/unpacked_useq/GRO11735_GRCh38_cellranger"

conda_sh="/hpc/dla_lti/yermanos_group/software/venvs4/miniforge3/etc/profile.d/conda.sh"
r_env="/hpc/dla_lti/yermanos_group/software/venvs4/r"

cfg="${project_root}/00_config"
inputs="${project_root}/01_inputs"
results="${project_root}/02_results/vdj_build"
logs="${project_root}/03_logs"
scripts="${project_root}/04_scripts"

mkdir -p "${cfg}" "${inputs}" "${results}" "${logs}" "${scripts}"

job_tmp_id="${SLURM_JOB_ID:-$$}"
export TMPDIR="${TMPDIR:-/scratch/${job_tmp_id}}"
mkdir -p "${TMPDIR}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

echo "Job started: $(date)"
echo "HOST: $(hostname)"
echo "PROJECT: ${project_root}"
echo "GRO11679_ROOT: ${gro11679_root}"
echo "GRO11735_ROOT: ${gro11735_root}"
echo "TMPDIR: ${TMPDIR}"
echo "CPUS: ${SLURM_CPUS_PER_TASK}"

source "${conda_sh}"
conda activate "${r_env}"

r_script="${scripts}/01_vdj_build_per_sample_useq_trimmed.R"

cat > "${r_script}" <<'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Platypus)
})

project_root <- "/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-04_AF_Pipeline_with_useq_aligned_data"

source_roots <- c(
  "/hpc/dla_lti/yermanos_group/GRO11679/aligned_and_demultiplexed_by_useq_GRO11679/unpacked_useq/GRO11679_GRCh38_cellranger",
  "/hpc/dla_lti/yermanos_group/GRO11735/aligned_and_demultiplexed_by_useq_GRO11735/unpacked_useq/GRO11735_GRCh38_cellranger"
)

cfg <- file.path(project_root, "00_config")
results_root <- file.path(project_root, "02_results", "vdj_build")
samples_tsv <- file.path(cfg, "samples_per_sample_vdj_dirs.tsv")

dir.create(cfg, showWarnings = FALSE, recursive = TRUE)
dir.create(results_root, showWarnings = FALSE, recursive = TRUE)

must_have <- c(
  "filtered_contig_annotations.csv",
  "filtered_contig.fasta",
  "consensus_annotations.csv",
  "consensus.fasta",
  "clonotypes.csv",
  "concat_ref.fasta"
)

read_tsv_base <- function(path) {
  x <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (nrow(x) < 1) stop("Empty TSV: ", path)
  x
}

find_sample_vdj_dirs <- function(roots) {
  out <- character(0)
  for (root in roots) {
    if (!dir.exists(root)) {
      warning("Source root does not exist and will be skipped: ", root)
      next
    }
    hits <- Sys.glob(file.path(root, "*", "outs", "per_sample_outs", "*", "vdj_b"))
    hits <- hits[file.info(hits)$isdir %in% TRUE]
    out <- c(out, hits)
  }
  out <- unique(normalizePath(out, winslash = "/", mustWork = FALSE))
  out <- sort(out)
  out
}

assert_sample_vdj_dir <- function(d, sample_id) {
  if (!dir.exists(d)) stop("VDJ dir not found for ", sample_id, ": ", d)
  missing <- must_have[!file.exists(file.path(d, must_have))]
  if (length(missing) > 0) {
    stop(
      "VDJ dir for ", sample_id, " is missing files: ",
      paste(missing, collapse = ", "),
      "\nDir: ", d
    )
  }
}

parse_sample_metadata <- function(vdj_dir) {
  sample_folder <- basename(dirname(vdj_dir))
  run_dir <- basename(dirname(dirname(dirname(dirname(vdj_dir)))))

  source_project <- if (grepl("/GRO11679/", vdj_dir)) {
    "GRO11679"
  } else if (grepl("/GRO11735/", vdj_dir)) {
    "GRO11735"
  } else {
    "UNKNOWN"
  }

  patient_id <- NA_character_
  sample_id  <- NA_character_
  tissue     <- NA_character_
  side       <- NA_character_

  if (grepl("_Sc-Seq1_PBMC_Bcells_", sample_folder, fixed = TRUE)) {
    patient_id <- sub("_Sc-Seq1_PBMC_Bcells_.*$", "", sample_folder)
    sample_id  <- paste0(patient_id, "_B")
    tissue     <- "blood"
    side       <- "NA"
  } else if (grepl("-Sc-Seq1-", sample_folder, fixed = TRUE) &&
             grepl("-SFMC-Bcells-", sample_folder, fixed = TRUE)) {
    m <- regexec("^((HC|P)[0-9]+)-Sc-Seq1-([A-Z]+)-SFMC-Bcells-.*$", sample_folder)
    reg <- regmatches(sample_folder, m)[[1]]
    if (length(reg) < 4) {
      stop("Could not parse SF sample folder name: ", sample_folder)
    }
    patient_id <- reg[2]
    side_raw <- reg[4]

    if (identical(side_raw, "K")) {
      side <- "LK"
    } else {
      side <- side_raw
    }

    sample_id <- paste0(patient_id, "_", side)
    tissue <- "SF"
  } else {
    stop(
      "Could not parse sample metadata from folder name:\n",
      sample_folder, "\n",
      "Path: ", vdj_dir
    )
  }

  data.frame(
    source_project = source_project,
    run_dir = run_dir,
    sample_folder = sample_folder,
    patient_id = patient_id,
    sample_id = sample_id,
    tissue = tissue,
    side = side,
    vdj_dir = vdj_dir,
    stringsAsFactors = FALSE
  )
}

make_vdj_build_compat <- function() {
  fn <- get("VDJ_build", envir = asNamespace("Platypus"))

  biostrings_exports <- getNamespaceExports("Biostrings")

  if ("pairwiseAlignment" %in% biostrings_exports) {
    message("Biostrings still exports pairwiseAlignment. Using Platypus::VDJ_build unchanged.")
    return(fn)
  }

  if (!requireNamespace("pwalign", quietly = TRUE)) {
    stop(
      paste0(
        "Biostrings does not export pairwiseAlignment in this environment, ",
        "and pwalign is not installed.\n",
        "Install pwalign in the active R environment and rerun.\n"
      )
    )
  }

  fn_txt <- paste(deparse(fn), collapse = "\n")

  if (!grepl("Biostrings::pairwiseAlignment", fn_txt, fixed = TRUE)) {
    stop(
      "Could not find 'Biostrings::pairwiseAlignment' inside Platypus::VDJ_build. ",
      "The runtime compatibility patch cannot be applied safely."
    )
  }

  fn_txt_patched <- gsub(
    "Biostrings::pairwiseAlignment",
    "pwalign::pairwiseAlignment",
    fn_txt,
    fixed = TRUE
  )

  fn2 <- eval(parse(text = fn_txt_patched))
  environment(fn2) <- environment(fn)

  message("Using patched VDJ_build: Biostrings::pairwiseAlignment -> pwalign::pairwiseAlignment")
  fn2
}

VDJ_build_compat <- make_vdj_build_compat()

if (!file.exists(samples_tsv) || file.info(samples_tsv)$size == 0) {
  vdj_dirs <- find_sample_vdj_dirs(source_roots)
  if (length(vdj_dirs) < 1) {
    stop("No per sample vdj_b folders were found under the provided source roots.")
  }

  smp_list <- lapply(vdj_dirs, parse_sample_metadata)
  smp <- do.call(rbind, smp_list)

  smp <- smp[order(smp$sample_id, smp$source_project, smp$run_dir), , drop = FALSE]

  dup <- duplicated(smp$sample_id)
  if (any(dup)) {
    stop(
      "Duplicate sample_id values were produced while auto building the config: ",
      paste(unique(smp$sample_id[dup]), collapse = ", "),
      "\nPlease edit ", samples_tsv, " manually."
    )
  }

  write.table(
    smp,
    file = samples_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("Created sample config:", samples_tsv, "\n")
} else {
  cat("Using existing sample config:", samples_tsv, "\n")
}

smp <- read_tsv_base(samples_tsv)

required_cols <- c(
  "source_project",
  "run_dir",
  "sample_folder",
  "patient_id",
  "sample_id",
  "tissue",
  "side",
  "vdj_dir"
)

missing_cols <- setdiff(required_cols, colnames(smp))
if (length(missing_cols) > 0) {
  stop(
    "Sample config is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\nFile: ", samples_tsv
  )
}

if (anyDuplicated(smp$sample_id)) {
  stop(
    "sample_id values are not unique in: ", samples_tsv,
    "\nDuplicates: ",
    paste(unique(smp$sample_id[duplicated(smp$sample_id)]), collapse = ", ")
  )
}

cat("R version:", R.version.string, "\n")
cat("Platypus version:", as.character(utils::packageVersion("Platypus")), "\n")
cat("Biostrings version:", as.character(utils::packageVersion("Biostrings")), "\n")
if (requireNamespace("pwalign", quietly = TRUE)) {
  cat("pwalign version:", as.character(utils::packageVersion("pwalign")), "\n")
} else {
  cat("pwalign not installed in this environment\n")
}

qc_lines <- list()

for (i in seq_len(nrow(smp))) {
  source_project <- smp$source_project[i]
  run_dir        <- smp$run_dir[i]
  sample_folder  <- smp$sample_folder[i]
  patient_id     <- smp$patient_id[i]
  sample_id      <- smp$sample_id[i]
  tissue         <- smp$tissue[i]
  side           <- smp$side[i]
  vdj_dir        <- smp$vdj_dir[i]

  cat("\n============================================================\n")
  cat("Running VDJ_build for sample:", sample_id, "\n")
  cat("Source project:", source_project, "\n")
  cat("Run dir:", run_dir, "\n")
  cat("Sample folder:", sample_folder, "\n")
  cat("VDJ dir:", vdj_dir, "\n")

  assert_sample_vdj_dir(vdj_dir, sample_id)

  vdj_s <- VDJ_build_compat(
    VDJ.sample.list        = c(vdj_dir),
    remove.divergent.cells = TRUE,
    complete.cells.only    = TRUE,
    trim.germlines         = TRUE
  )

  vdj_s$source_project <- source_project
  vdj_s$run_dir <- run_dir
  vdj_s$sample_folder <- sample_folder
  vdj_s$patient_id <- patient_id
  vdj_s$sample_id <- sample_id
  vdj_s$tissue <- tissue
  vdj_s$side <- side

  if ("clonotype_id" %in% colnames(vdj_s)) {
    vdj_s$clonotype_id_10x_raw <- vdj_s$clonotype_id
    vdj_s$clonotype_id <- paste0(sample_id, "__", vdj_s$clonotype_id)
  }

  out_dir <- file.path(results_root, "by_sample", sample_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  out_rds <- file.path(out_dir, paste0("VDJ_build_", sample_id, ".rds.xz"))
  saveRDS(vdj_s, out_rds, compress = "xz")

  n_barcodes <- if ("barcode" %in% colnames(vdj_s)) {
    length(unique(vdj_s$barcode))
  } else {
    NA_integer_
  }

  n_clonotypes <- if ("clonotype_id" %in% colnames(vdj_s)) {
    length(unique(vdj_s$clonotype_id))
  } else {
    NA_integer_
  }

  qc_lines[[length(qc_lines) + 1]] <- data.frame(
    source_project = source_project,
    run_dir = run_dir,
    sample_folder = sample_folder,
    patient_id = patient_id,
    sample_id = sample_id,
    tissue = tissue,
    side = side,
    n_rows = nrow(vdj_s),
    n_unique_barcodes = n_barcodes,
    n_unique_clonotypes = n_clonotypes,
    vdj_dir = vdj_dir,
    stringsAsFactors = FALSE
  )

  cat("Saved sample:", sample_id, "rows:", nrow(vdj_s), "\n")
}

qc <- do.call(rbind, qc_lines)
qc_path <- file.path(results_root, "QC_counts_by_sample.tsv")
write.table(qc, qc_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWrote QC table:", qc_path, "\n")

patients <- unique(smp$patient_id)

for (pid in patients) {
  smp_ids <- smp$sample_id[smp$patient_id == pid]
  smp_ids <- unique(smp_ids)

  paths <- file.path(
    results_root,
    "by_sample",
    smp_ids,
    paste0("VDJ_build_", smp_ids, ".rds.xz")
  )

  ok <- file.exists(paths)
  if (any(!ok)) {
    stop(
      "Missing per sample RDS for patient ", pid, ": ",
      paste(smp_ids[!ok], collapse = ", ")
    )
  }

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
