#!/bin/bash
#SBATCH --job-name=enclone_reclonotype_by_patient
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-04_AF_Pipeline_with_useq_aligned_data/03_logs/%x_%j.out
#SBATCH --error=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-04_AF_Pipeline_with_useq_aligned_data/03_logs/%x_%j.err

set -euo pipefail

project_root="/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-04_AF_Pipeline_with_useq_aligned_data"

conda_sh="/hpc/dla_lti/yermanos_group/software/venvs4/miniforge3/etc/profile.d/conda.sh"
r_env="/hpc/dla_lti/yermanos_group/software/venvs4/r"
enclone_bin="/home/dla_lti/egrossewichtrup/bin/enclone"

cfg="${project_root}/00_config"
results="${project_root}/02_results/enclone_reclonotype"
logs="${project_root}/03_logs"
scripts="${project_root}/04_scripts"

mkdir -p "${cfg}" "${results}" "${logs}" "${scripts}"

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
echo "TMPDIR: ${TMPDIR}"
echo "CPUS: ${SLURM_CPUS_PER_TASK}"
echo "ENCLONE: ${enclone_bin}"

if [[ ! -x "${enclone_bin}" ]]; then
  echo "ERROR: enclone binary not found or not executable: ${enclone_bin}" >&2
  exit 1
fi

source "${conda_sh}"
conda activate "${r_env}"

export PATH="$(dirname "${enclone_bin}"):${PATH}"

echo "R path: $(which R)"
echo "R version:"
R --version | head -n 1

echo "enclone version:"
"${enclone_bin}" version

r_script="${scripts}/02_enclone_reclonotype_by_patient.R"

cat > "${r_script}" <<'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Platypus)
  library(dplyr)
  library(stringr)
})

project_root <- "/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-04_AF_Pipeline_with_useq_aligned_data"
enclone_bin  <- "/home/dla_lti/egrossewichtrup/bin/enclone"

cfg              <- file.path(project_root, "00_config")
vdj_build_root   <- file.path(project_root, "02_results", "vdj_build")
enclone_root     <- file.path(project_root, "02_results", "enclone_reclonotype")
reorg_root       <- file.path(enclone_root, "cellranger_per_sample")
by_patient_root  <- file.path(enclone_root, "by_patient")
by_sample_root   <- file.path(enclone_root, "by_sample_from_patient")
manifest_root    <- file.path(enclone_root, "manifests")
debug_root       <- file.path(enclone_root, "debug_enclone_files")

dir.create(enclone_root, recursive = TRUE, showWarnings = FALSE)
dir.create(reorg_root, recursive = TRUE, showWarnings = FALSE)
dir.create(by_patient_root, recursive = TRUE, showWarnings = FALSE)
dir.create(by_sample_root, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_root, recursive = TRUE, showWarnings = FALSE)
dir.create(debug_root, recursive = TRUE, showWarnings = FALSE)

samples_tsv <- file.path(cfg, "samples_per_sample_vdj_dirs.tsv")
if (!file.exists(samples_tsv)) {
  stop("Missing config file: ", samples_tsv)
}

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

samples <- read_tsv_base(samples_tsv)

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

missing_cols <- setdiff(required_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop(
    "Sample config is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\nFile: ", samples_tsv
  )
}

if (anyDuplicated(samples$sample_id)) {
  stop(
    "sample_id values are not unique in ",
    samples_tsv,
    "\nDuplicates: ",
    paste(unique(samples$sample_id[duplicated(samples$sample_id)]), collapse = ", ")
  )
}

must_have_sample_files <- c(
  "filtered_contig_annotations.csv",
  "filtered_contig.fasta",
  "consensus_annotations.csv",
  "consensus.fasta",
  "clonotypes.csv",
  "concat_ref.fasta"
)

must_have_enclone_files <- c(
  must_have_sample_files,
  "all_contig_annotations.json"
)

normalize_path <- function(x) {
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

run_dir_from_vdj_dir <- function(vdj_dir) {
  normalize_path(dirname(dirname(dirname(dirname(vdj_dir)))))
}

find_all_contig_json <- function(vdj_dir) {
  run_dir <- run_dir_from_vdj_dir(vdj_dir)

  candidate_paths <- c(
    file.path(vdj_dir, "all_contig_annotations.json"),
    file.path(run_dir, "outs", "per_sample_outs", basename(dirname(vdj_dir)), "vdj_b", "all_contig_annotations.json"),
    file.path(run_dir, "outs", "multi", "vdj_b", "all_contig_annotations.json"),
    file.path(run_dir, "outs", "vdj_b", "all_contig_annotations.json")
  )

  candidate_paths <- unique(candidate_paths[file.exists(candidate_paths)])

  if (length(candidate_paths) > 0) {
    return(normalize_path(candidate_paths[1]))
  }

  recursive_hits <- list.files(
    path = run_dir,
    pattern = "^all_contig_annotations\\.json$",
    recursive = TRUE,
    full.names = TRUE
  )

  recursive_hits <- recursive_hits[file.exists(recursive_hits)]
  recursive_hits <- unique(normalize_path(recursive_hits))

  if (length(recursive_hits) > 0) {
    preferred_multi <- recursive_hits[grepl("/outs/multi/vdj_b/all_contig_annotations\\.json$", recursive_hits)]
    if (length(preferred_multi) > 0) {
      return(preferred_multi[1])
    }
    return(recursive_hits[1])
  }

  stop(
    "Could not find all_contig_annotations.json for vdj_dir:\n",
    vdj_dir, "\n",
    "Searched under run dir:\n",
    run_dir
  )
}

copy_required_file <- function(src, dest) {
  ok <- file.copy(
    from = src,
    to = dest,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  if (!isTRUE(ok) || !file.exists(dest)) {
    stop("Could not copy file:\n  from: ", src, "\n  to:   ", dest)
  }
}

stage_sample_for_enclone <- function(sample_id, vdj_dir) {
  if (!dir.exists(vdj_dir)) {
    stop("vdj_dir does not exist for sample ", sample_id, ": ", vdj_dir)
  }

  dest_dir <- file.path(reorg_root, sample_id)

  if (dir.exists(dest_dir)) {
    unlink(dest_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  missing_sample_files <- must_have_sample_files[!file.exists(file.path(vdj_dir, must_have_sample_files))]
  if (length(missing_sample_files) > 0) {
    stop(
      "Sample ", sample_id, " is missing required per sample files in ", vdj_dir, ":\n",
      paste(missing_sample_files, collapse = ", ")
    )
  }

  for (fname in must_have_sample_files) {
    copy_required_file(
      src = file.path(vdj_dir, fname),
      dest = file.path(dest_dir, fname)
    )
  }

  json_src <- find_all_contig_json(vdj_dir)
  copy_required_file(
    src = json_src,
    dest = file.path(dest_dir, "all_contig_annotations.json")
  )

  missing_after_stage <- must_have_enclone_files[!file.exists(file.path(dest_dir, must_have_enclone_files))]
  if (length(missing_after_stage) > 0) {
    stop(
      "Staging failed for sample ", sample_id, ". Missing in staged dir ", dest_dir, ":\n",
      paste(missing_after_stage, collapse = ", ")
    )
  }

  data.frame(
    sample_id = sample_id,
    source_vdj_dir = normalize_path(vdj_dir),
    staged_dir = normalize_path(dest_dir),
    json_source = json_src,
    stringsAsFactors = FALSE
  )
}

cat("Rebuilding Cell Ranger VDJ folders for enclone\n")
stage_manifest <- do.call(
  rbind,
  lapply(seq_len(nrow(samples)), function(i) {
    stage_sample_for_enclone(
      sample_id = samples$sample_id[i],
      vdj_dir = samples$vdj_dir[i]
    )
  })
)

write.table(
  stage_manifest,
  file = file.path(manifest_root, "sample_stage_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

normalize_barcode <- function(x) {
  x <- as.character(x)
  x <- sub("^.*_", "", x)
  x <- sub("-[0-9]+$", "", x)
  x
}

build_sample_groups <- function(vdj_df, samples.to.combine = NULL) {
  repertoire.number <- unique(as.character(vdj_df$sample_id))
  groups <- list()

  if (is.null(samples.to.combine)) {
    for (sid in repertoire.number) {
      groups[[length(groups) + 1]] <- vdj_df[vdj_df$sample_id == sid, , drop = FALSE]
    }
    return(groups)
  }

  used <- character(0)

  if (is.list(samples.to.combine)) {
    for (grp in samples.to.combine) {
      grp <- unique(as.character(grp))
      if (!all(grp %in% repertoire.number)) {
        stop(
          "Some samples in samples.to.combine were not found in VDJ: ",
          paste(setdiff(grp, repertoire.number), collapse = ", ")
        )
      }
      groups[[length(groups) + 1]] <- vdj_df[vdj_df$sample_id %in% grp, , drop = FALSE]
      used <- c(used, grp)
    }
  } else {
    grp <- unique(as.character(samples.to.combine))
    if (!all(grp %in% repertoire.number)) {
      stop(
        "Some samples in samples.to.combine were not found in VDJ: ",
        paste(setdiff(grp, repertoire.number), collapse = ", ")
      )
    }
    groups[[length(groups) + 1]] <- vdj_df[vdj_df$sample_id %in% grp, , drop = FALSE]
    used <- c(used, grp)
  }

  remaining <- setdiff(repertoire.number, used)
  for (sid in remaining) {
    groups[[length(groups) + 1]] <- vdj_df[vdj_df$sample_id == sid, , drop = FALSE]
  }

  groups
}

write_meta_csv_unquoted <- function(path, sample_names, vdj_directory) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("bcr,origin", con)

  lines <- vapply(
    sample_names,
    function(sid) paste(file.path(vdj_directory, sid), sid, sep = ","),
    character(1)
  )

  writeLines(lines, con)
}

VDJ_enclone_demultiplexed_safe <- function(
  VDJ,
  VDJ.directory,
  samples.to.clonotype = "all",
  samples.to.combine = NULL,
  enclone.binary,
  debug.dir
) {
  if (missing(VDJ) || is.null(VDJ)) {
    stop("Please provide input data as VDJ")
  }
  if (missing(VDJ.directory) || is.null(VDJ.directory)) {
    stop("Please provide VDJ.directory")
  }
  if (missing(enclone.binary) || !file.exists(enclone.binary)) {
    stop("enclone binary not found: ", enclone.binary)
  }
  if (missing(debug.dir) || is.null(debug.dir)) {
    stop("Please provide debug.dir")
  }
  if (!("sample_id" %in% colnames(VDJ))) {
    stop("VDJ must contain a sample_id column")
  }
  if (!("barcode" %in% colnames(VDJ))) {
    stop("VDJ must contain a barcode column")
  }
  if (!("celltype" %in% colnames(VDJ))) {
    VDJ$celltype <- "B cell"
  }
  if (!("X" %in% colnames(VDJ))) {
    VDJ$X <- seq_len(nrow(VDJ))
  }

  dir.create(debug.dir, recursive = TRUE, showWarnings = FALSE)

  repertoire.number <- unique(as.character(VDJ$sample_id))

  if (length(samples.to.clonotype) == 1 && identical(samples.to.clonotype, "all")) {
    samples_to_clonotype <- repertoire.number
  } else {
    samples_to_clonotype <- unique(as.character(samples.to.clonotype))
  }

  sample_dirs <- list.files(VDJ.directory)
  missing_sample_dirs <- setdiff(samples_to_clonotype, sample_dirs)
  if (length(missing_sample_dirs) > 0) {
    stop(
      "Unable to find staged VDJ folders for: ",
      paste(missing_sample_dirs, collapse = ", "),
      "\nVDJ.directory = ", VDJ.directory
    )
  }

  samples_not_clonotyped <- VDJ[!(VDJ$sample_id %in% samples_to_clonotype), , drop = FALSE]
  VDJ_work <- VDJ[VDJ$sample_id %in% samples_to_clonotype, , drop = FALSE]

  if (nrow(samples_not_clonotyped) > 0 && !("clonotype_id_10x" %in% colnames(samples_not_clonotyped))) {
    samples_not_clonotyped$clonotype_id_10x <- samples_not_clonotyped$clonotype_id
  }

  get_enclone_output <- function(sample_df) {
    sample_names <- unique(as.character(sample_df$sample_id))
    message("Processing sample group: ", paste(sample_names, collapse = ", "))

    file_stub <- paste(sample_names, collapse = "_")
    temp_file <- file.path(debug.dir, paste0(file_stub, ".csv"))
    meta_file <- file.path(debug.dir, paste0(file_stub, "_META.csv"))
    log_file  <- file.path(debug.dir, paste0(file_stub, "_enclone_stdout_stderr.txt"))

    celltype <- unique(as.character(sample_df$celltype))
    celltype <- celltype[!is.na(celltype)]
    if (length(celltype) != 1) {
      stop(
        "enclone only works on a single celltype per group. Found: ",
        paste(celltype, collapse = ", ")
      )
    }

    if (!identical(celltype, "B cell")) {
      stop("This script is currently set up for B cell data only. Found celltype: ", celltype)
    }

    for (sid in sample_names) {
      needed_here <- c("all_contig_annotations.json", "filtered_contig_annotations.csv", "clonotypes.csv", "concat_ref.fasta")
      miss <- needed_here[!file.exists(file.path(VDJ.directory, sid, needed_here))]
      if (length(miss) > 0) {
        stop(
          "Before enclone, staged sample folder is missing files for ", sid, ":\n",
          paste(miss, collapse = ", "), "\n",
          "Folder: ", file.path(VDJ.directory, sid)
        )
      }
    }

    write_meta_csv_unquoted(
      path = meta_file,
      sample_names = sample_names,
      vdj_directory = VDJ.directory
    )

    args <- c(
      paste0("META=", meta_file),
      paste0("POUT=", temp_file),
      "PCELL",
      "PCOLS=barcode,datasets,group_id",
      "NOPRINT",
      "ACCEPT_REUSE"
    )

    out <- tryCatch(
      system2(enclone.binary, args = args, stdout = TRUE, stderr = TRUE),
      error = function(e) e
    )

    if (!inherits(out, "error")) {
      writeLines(out, log_file)
    }

    if (inherits(out, "error")) {
      stop(
        "enclone failed for sample group ",
        paste(sample_names, collapse = ", "),
        "\n", conditionMessage(out),
        "\nMETA file kept at: ", meta_file
      )
    }

    status <- attr(out, "status")
    if (!is.null(status) && status != 0) {
      stop(
        "enclone exited with status ", status,
        " for sample group ", paste(sample_names, collapse = ", "),
        "\n",
        paste(out, collapse = "\n"),
        "\nMETA file kept at: ", meta_file,
        "\nPOUT target: ", temp_file,
        "\nLog file: ", log_file
      )
    }

    if (!file.exists(temp_file) || file.info(temp_file)$size <= 0) {
      stop(
        "enclone did not create parseable output for sample group ",
        paste(sample_names, collapse = ", "),
        "\nMETA file kept at: ", meta_file,
        "\nLog file: ", log_file
      )
    }

    list(
      sample_df = sample_df,
      sample_names = sample_names,
      csv_file = temp_file,
      meta_file = meta_file
    )
  }

  merge_enclone_clonotypes <- function(enclone_out) {
    new_clonotypes <- utils::read.csv(
      enclone_out$csv_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    needed <- c("barcode", "group_id", "datasets")
    if (!all(needed %in% colnames(new_clonotypes))) {
      stop(
        "Parseable enclone output is missing required columns in ",
        enclone_out$csv_file,
        "\nFound columns: ",
        paste(colnames(new_clonotypes), collapse = ", ")
      )
    }

    names(new_clonotypes)[names(new_clonotypes) == "group_id"] <- "new_group_id"
    names(new_clonotypes)[names(new_clonotypes) == "barcode"] <- "barcode_enclone"
    names(new_clonotypes)[names(new_clonotypes) == "datasets"] <- "sample_id_enclone"

    new_clonotypes$sample_id_enclone <- as.character(new_clonotypes$sample_id_enclone)
    new_clonotypes$sample_id_enclone <- sub(",.*$", "", new_clonotypes$sample_id_enclone)
    new_clonotypes$stripped_barcodes <- normalize_barcode(new_clonotypes$barcode_enclone)

    new_clonotypes <- new_clonotypes %>%
      dplyr::group_by(sample_id_enclone, stripped_barcodes) %>%
      dplyr::summarise(
        new_group_id = paste(new_group_id, collapse = "_"),
        .groups = "drop"
      )

    sample_df <- enclone_out$sample_df
    sample_df$stripped_barcodes <- normalize_barcode(sample_df$barcode)

    sample_out <- dplyr::left_join(
      sample_df,
      new_clonotypes,
      by = c("sample_id" = "sample_id_enclone", "stripped_barcodes" = "stripped_barcodes")
    )

    if ("clonotype_id_10x" %in% colnames(sample_out)) {
      sample_out$clonotype_id_before_enclone <- sample_out$clonotype_id
    } else {
      sample_out$clonotype_id_10x <- sample_out$clonotype_id
    }

    sample_out$clonotype_id <- ifelse(
      is.na(sample_out$new_group_id),
      "clonotypeNA",
      paste0(
        "clonotype",
        vapply(
          strsplit(sample_out$new_group_id, "_", fixed = TRUE),
          `[`,
          character(1),
          1L
        )
      )
    )

    freq_tab <- table(sample_out$clonotype_id, useNA = "no")
    sample_out$clonotype_frequency <- as.integer(freq_tab[sample_out$clonotype_id])

    sample_out$stripped_barcodes <- NULL
    sample_out$new_group_id <- NULL

    clonotype_map <- sample_out %>%
      dplyr::filter(clonotype_id != "clonotypeNA") %>%
      dplyr::distinct(clonotype_id, clonotype_frequency) %>%
      dplyr::arrange(dplyr::desc(clonotype_frequency), clonotype_id) %>%
      dplyr::mutate(new_clonotype_id = paste0("clonotype", dplyr::row_number()))

    sample_out <- sample_out %>%
      dplyr::left_join(
        clonotype_map,
        by = c("clonotype_id", "clonotype_frequency")
      ) %>%
      dplyr::mutate(
        clonotype_id = ifelse(
          clonotype_id == "clonotypeNA",
          "clonotypeNA",
          new_clonotype_id
        )
      )

    sample_out$new_clonotype_id <- NULL
    sample_out <- sample_out[!duplicated(sample_out$X), , drop = FALSE]
    sample_out <- sample_out[order(sample_out$X), , drop = FALSE]

    if (any(sample_out$clonotype_id == "clonotypeNA")) {
      message("NA clonotypes found, these remain labeled as clonotypeNA.")
    }

    sample_out
  }

  sample_dfs <- build_sample_groups(VDJ_work, samples.to.combine = samples.to.combine)

  enclone_outs <- lapply(sample_dfs, get_enclone_output)
  sample_outs  <- lapply(enclone_outs, merge_enclone_clonotypes)

  VDJ_out <- dplyr::bind_rows(sample_outs, samples_not_clonotyped)
  VDJ_out <- VDJ_out[order(VDJ_out$X), , drop = FALSE]
  VDJ_out$X <- NULL

  VDJ_out
}

patient_ids <- unique(as.character(samples$patient_id))

patient_qc <- list()
sample_qc  <- list()
all_patient_results <- list()

cat("R version: ", R.version.string, "\n", sep = "")
cat("Platypus version: ", as.character(utils::packageVersion("Platypus")), "\n", sep = "")
cat("dplyr version: ", as.character(utils::packageVersion("dplyr")), "\n", sep = "")
cat("stringr version: ", as.character(utils::packageVersion("stringr")), "\n", sep = "")

for (pid in patient_ids) {
  patient_file <- file.path(
    vdj_build_root,
    "by_patient",
    pid,
    paste0("VDJ_build_patient_", pid, ".rds.xz")
  )

  if (!file.exists(patient_file)) {
    stop("Missing patient RDS from step 1: ", patient_file)
  }

  vdj_patient <- readRDS(patient_file)

  if (!("sample_id" %in% colnames(vdj_patient))) {
    stop("Patient file lacks sample_id column: ", patient_file)
  }
  if (!("barcode" %in% colnames(vdj_patient))) {
    stop("Patient file lacks barcode column: ", patient_file)
  }
  if (!("patient_id" %in% colnames(vdj_patient))) {
    vdj_patient$patient_id <- pid
  }
  if (!("celltype" %in% colnames(vdj_patient))) {
    vdj_patient$celltype <- "B cell"
  }

  patient_sample_ids <- unique(as.character(samples$sample_id[samples$patient_id == pid]))
  patient_sample_ids <- patient_sample_ids[patient_sample_ids %in% unique(as.character(vdj_patient$sample_id))]

  if (length(patient_sample_ids) < 1) {
    stop("No sample_ids found for patient ", pid, " in ", patient_file)
  }

  cat("\n============================================================\n")
  cat("Running enclone for patient: ", pid, "\n", sep = "")
  cat("Samples: ", paste(patient_sample_ids, collapse = ", "), "\n", sep = "")
  cat("Input file: ", patient_file, "\n", sep = "")

  combine_arg <- if (length(patient_sample_ids) > 1) patient_sample_ids else NULL
  patient_debug_dir <- file.path(debug_root, pid)

  vdj_enclone <- VDJ_enclone_demultiplexed_safe(
    VDJ = vdj_patient,
    VDJ.directory = reorg_root,
    samples.to.clonotype = patient_sample_ids,
    samples.to.combine = combine_arg,
    enclone.binary = enclone_bin,
    debug.dir = patient_debug_dir
  )

  vdj_enclone$patient_id <- pid
  vdj_enclone$clonotype_id_within_patient <- vdj_enclone$clonotype_id
  vdj_enclone$clonotype_id <- ifelse(
    vdj_enclone$clonotype_id == "clonotypeNA",
    "clonotypeNA",
    paste0(pid, "__", vdj_enclone$clonotype_id_within_patient)
  )

  out_dir_patient <- file.path(by_patient_root, pid)
  dir.create(out_dir_patient, recursive = TRUE, showWarnings = FALSE)

  patient_out_rds <- file.path(
    out_dir_patient,
    paste0("VDJ_enclone_patient_", pid, ".rds.xz")
  )

  saveRDS(vdj_enclone, patient_out_rds, compress = "xz")

  n_patient_barcodes <- length(unique(vdj_enclone$barcode))
  n_patient_clonotypes <- length(unique(vdj_enclone$clonotype_id[vdj_enclone$clonotype_id != "clonotypeNA"]))

  patient_qc[[length(patient_qc) + 1]] <- data.frame(
    patient_id = pid,
    sample_ids = paste(patient_sample_ids, collapse = ","),
    n_samples = length(patient_sample_ids),
    n_rows = nrow(vdj_enclone),
    n_unique_barcodes = n_patient_barcodes,
    n_unique_clonotypes = n_patient_clonotypes,
    output_rds = patient_out_rds,
    stringsAsFactors = FALSE
  )

  for (sid in patient_sample_ids) {
    vdj_sample <- vdj_enclone[vdj_enclone$sample_id == sid, , drop = FALSE]

    out_dir_sample <- file.path(by_sample_root, sid)
    dir.create(out_dir_sample, recursive = TRUE, showWarnings = FALSE)

    sample_out_rds <- file.path(
      out_dir_sample,
      paste0("VDJ_enclone_", sid, ".rds.xz")
    )

    saveRDS(vdj_sample, sample_out_rds, compress = "xz")

    n_sample_barcodes <- length(unique(vdj_sample$barcode))
    n_sample_clonotypes <- length(unique(vdj_sample$clonotype_id[vdj_sample$clonotype_id != "clonotypeNA"]))

    sample_qc[[length(sample_qc) + 1]] <- data.frame(
      patient_id = pid,
      sample_id = sid,
      n_rows = nrow(vdj_sample),
      n_unique_barcodes = n_sample_barcodes,
      n_unique_clonotypes = n_sample_clonotypes,
      output_rds = sample_out_rds,
      stringsAsFactors = FALSE
    )
  }

  all_patient_results[[pid]] <- vdj_enclone

  cat("Saved patient output: ", patient_out_rds, "\n", sep = "")
  cat("Rows: ", nrow(vdj_enclone), "\n", sep = "")
  cat("Unique barcodes: ", n_patient_barcodes, "\n", sep = "")
  cat("Unique clonotypes: ", n_patient_clonotypes, "\n", sep = "")
}

all_patients <- dplyr::bind_rows(all_patient_results)

all_patients_rds <- file.path(enclone_root, "all_patients_reclonotyped.rds.xz")
saveRDS(all_patients, all_patients_rds, compress = "xz")

patient_qc_df <- dplyr::bind_rows(patient_qc)
sample_qc_df  <- dplyr::bind_rows(sample_qc)

write.table(
  patient_qc_df,
  file = file.path(manifest_root, "QC_counts_by_patient.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_qc_df,
  file = file.path(manifest_root, "QC_counts_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nDONE\n")
cat("Combined output: ", all_patients_rds, "\n", sep = "")
cat("Patient QC: ", file.path(manifest_root, "QC_counts_by_patient.tsv"), "\n", sep = "")
cat("Sample QC: ", file.path(manifest_root, "QC_counts_by_sample.tsv"), "\n", sep = "")
cat("Stage manifest: ", file.path(manifest_root, "sample_stage_manifest.tsv"), "\n", sep = "")
cat("Debug files: ", debug_root, "\n", sep = "")
RSCRIPT

chmod +x "${r_script}"

echo "Running R script: ${r_script}"
Rscript "${r_script}"

echo "Job finished: $(date)"