#!/usr/bin/env bash
#SBATCH --job-name=JIA_VDJ_ownCR
#SBATCH --partition=cpu
#SBATCH --account=dla_lti
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-10_VDJ_build_own_CR/03_logs/%x_%j.out
#SBATCH --error=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-10_VDJ_build_own_CR/03_logs/%x_%j.err

set -Eeuo pipefail
trap 'echo "ERROR: shell command failed (line ${LINENO}, exit ${?})." >&2' ERR

# JIA unified VDJ_build v1.3. Preserves v1.2 partial outputs.
# Create 03_logs before sbatch. This sbatch generates a separate versioned R file.
PROJECT_ROOT="/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-10_VDJ_build_own_CR"
MAP="${PROJECT_ROOT}/00_config/samples_JIA_16.tsv"
CONDA_SH="/hpc/dla_lti/yermanos_group/software/venvs4/miniforge3/etc/profile.d/conda.sh"
R_ENV="/hpc/dla_lti/yermanos_group/software/venvs4/r"
R_SCRIPT="${PROJECT_ROOT}/04_scripts/01_unified_jia_vdj_build_v1_3.R"
MODE="${1:---run}"
[[ "${MODE}" == "--run" || "${MODE}" == "--check" ]] || { echo 'Usage: bash script.sbatch [--check|--run]' >&2; exit 2; }
[[ -d "${PROJECT_ROOT}/03_logs" ]] || { echo "Create ${PROJECT_ROOT}/03_logs BEFORE sbatch." >&2; exit 2; }
[[ -s "${MAP}" ]] || { echo "Missing mapping TSV: ${MAP}" >&2; exit 2; }
[[ -f "${CONDA_SH}" ]] || { echo "Missing conda activation script: ${CONDA_SH}" >&2; exit 2; }
mkdir -p "${PROJECT_ROOT}/04_scripts" "${PROJECT_ROOT}/01_tmp"

# Never silently replace a generated R source from an earlier script version.
TEMP_R="$(mktemp "${PROJECT_ROOT}/04_scripts/.vdj_source_XXXXXXXX")"
trap 'rm -f "${TEMP_R}"' EXIT
cat > "${TEMP_R}" <<'RSCRIPT'
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[[1L]] %in% c("--check", "--run")) {
  stop("Usage: Rscript 01_unified_jia_vdj_build.R [--check|--run] PROJECT_ROOT MAPPING_TSV")
}
mode <- args[[1L]]
project_root <- args[[2L]]
map_path <- args[[3L]]
# Keep the failed v1.2 results and its R script untouched. A rerun of v1.3
# also refuses to overwrite its own results; archive before another attempt.
results_root <- file.path(project_root, "02_results", "vdj_build_v1_3")

say <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
fail <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)
valid <- function(x) !is.na(x) & nzchar(trimws(as.character(x))) & !(as.character(x) %in% c("None", "NULL", "NA"))
clean_missing <- function(x) { x <- as.character(x); x[!valid(x)] <- NA_character_; x }
write_tsv <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE,
                                              row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

# A frozen expected biological mapping; the TSV is editable for PATHS but cannot
# silently redefine the cohort. Changing cohort composition requires code review.
expected <- read.delim(text = paste(c(
  "source_project\tcapture_no\tsample_id\tpatient_id\ttissue\tside",
  "GRO11679\tC1\tHC001_B\tHC001\tblood\tNA",
  "GRO11679\tC1\tHC002_B\tHC002\tblood\tNA",
  "GRO11679\tC6\tHC003_B\tHC003\tblood\tNA",
  "GRO11679\tC6\tP003_B\tP003\tblood\tNA",
  "GRO11679\tC8\tP005_B\tP005\tblood\tNA",
  "GRO11679\tC8\tP006_B\tP006\tblood\tNA",
  "GRO11679\tC9\tP002_B\tP002\tblood\tNA",
  "GRO11679\tC9\tP007_B\tP007\tblood\tNA",
  "GRO11679\tC9\tP008_B\tP008\tblood\tNA",
  "GRO11679\tC10\tP009_B\tP009\tblood\tNA",
  "GRO11679\tC10\tP010_B\tP010\tblood\tNA",
  "GRO11735\tC2\tP001_LK\tP001\tSF\tLK",
  "GRO11735\tC3\tP001_RK\tP001\tSF\tRK",
  "GRO11735\tC4\tP002_LK\tP002\tSF\tLK",
  "GRO11735\tC5\tP002_RK\tP002\tSF\tRK",
  "GRO11735\tC7\tP004_LK\tP004\tSF\tLK"
), collapse = "\n"), sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, check.names = FALSE)
meta_cols <- c("source_project", "capture_no", "sample_id", "patient_id", "tissue", "side")
required_cols <- c(meta_cols, "run_dir", "vdj_dir")
needed_files <- c("filtered_contig_annotations.csv", "filtered_contig.fasta",
                  "consensus_annotations.csv", "consensus.fasta",
                  "clonotypes.csv", "concat_ref.fasta")

# Cell Ranger's CSV 'sample' can be a technical capture/run label even inside
# per_sample_outs/<biological_sample>/vdj_b. Do not confuse this column with
# the authoritative biological sample identity recorded in the validated TSV.
# Accept only the known biological labels OR the exact configured run name.
check_csv_sample_label <- function(x, s) {
  sid <- s$sample_id[[1L]]
  if (!("sample" %in% names(x))) return(NA_character_)
  labels <- unique(trimws(clean_missing(x$sample)))
  labels <- labels[!is.na(labels)]
  allowed <- unique(c(sid, s$patient_id[[1L]], basename(dirname(s$vdj_dir[[1L]])),
                      basename(s$run_dir[[1L]])))
  if (length(labels) > 1L || (length(labels) == 1L && !(labels[[1L]] %in% allowed))) {
    fail("[",sid,"] Cell Ranger CSV 'sample' is not a permitted biological/sample-folder/run label: ",
         paste(labels, collapse=", "), "; permitted: ",paste(allowed,collapse=", "))
  }
  if (length(labels) == 1L && identical(labels[[1L]], basename(s$run_dir[[1L]])))
    say("[",sid,"] Cell Ranger CSV 'sample' is technical run label: ",labels[[1L]],
        " (biological sample taken from validated mapping)")
  if (length(labels)) labels[[1L]] else NA_character_
}

preflight <- function() {
  if (!file.exists(map_path)) fail("Configuration not found: ", map_path)
  smp <- read.delim(map_path, sep = "\t", header = TRUE, quote = "\"",
                    comment.char = "", stringsAsFactors = FALSE, na.strings = c("NA", ""),
                    check.names = FALSE, colClasses = "character")
  if (!setequal(names(smp), required_cols) || anyDuplicated(names(smp)))
    fail("TSV columns must be exactly: ", paste(required_cols, collapse = ", "))
  if (nrow(smp) != 16L || anyDuplicated(smp$sample_id) || anyNA(smp$sample_id))
    fail("Expected exactly 16 nonmissing, distinct sample IDs; found ", nrow(smp))
  if (!setequal(smp$sample_id, expected$sample_id))
    fail("Sample set mismatch. Missing: ", paste(setdiff(expected$sample_id, smp$sample_id), collapse=","),
         "; extra: ", paste(setdiff(smp$sample_id, expected$sample_id), collapse=","))
  if (anyDuplicated(smp$vdj_dir) || anyNA(smp$vdj_dir) || any(!valid(smp$vdj_dir)))
    fail("Duplicated or missing vdj_dir in mapping")
  smp <- smp[match(expected$sample_id, smp$sample_id), , drop = FALSE]
  rownames(smp) <- NULL
  for (col in meta_cols) {
    got <- smp[[col]]; want <- as.character(expected[[col]])
    if (col == "side") want[is.na(expected[[col]])] <- NA_character_
    bad <- which(!(is.na(got) & is.na(want)) & (is.na(got) | is.na(want) | got != want))
    if (length(bad)) fail("Invalid ", col, " for ", paste(smp$sample_id[bad], collapse=", "),
                          "; expected mapping is fixed in the script")
  }
  if (length(unique(smp$patient_id)) != 13L || length(unique(smp$capture_no)) != 10L ||
      sum(smp$tissue == "blood") != 11L || sum(smp$tissue == "SF") != 5L)
    fail("Cohort counts inconsistent (13 patients / 10 captures / 11 blood / 5 SF)")
  # Within one hashed capture, a 10x droplet barcode must not appear in two
  # different per-sample VDJ tables. Detect accidental reuse of pooled input.
  seen_barcodes <- list()
  for (i in seq_len(nrow(smp))) {
    s <- smp[i, , drop=FALSE]; sid <- s$sample_id[[1L]]
    run <- s$run_dir[[1L]]; d <- s$vdj_dir[[1L]]
    if (!valid(run) || !valid(d) || !startsWith(run, "/") || !startsWith(d, "/"))
      fail("[", sid, "] Absolute, nonmissing run_dir and vdj_dir required")
    if (grepl("useq|unpacked|aligned_and_demultiplexed_by_useq", paste(run,d), ignore.case=TRUE))
      fail("[", sid, "] REFUSED: USEQ-derived directory: ", d)
    permitted <- if (s$source_project == "GRO11679")
      "/hpc/dla_lti/yermanos_group/GRO11679/2026-06_alignment_CR/cellranger_multi_outputs/"
      else "/hpc/dla_lti/yermanos_group/GRO11735/alignment_CR/"
    if (!startsWith(paste0(run, "/"), permitted) || !startsWith(paste0(d, "/"), permitted))
      fail("[", sid, "] Path outside allowlisted own Cell Ranger alignment root: ", d)
    if (!identical(dirname(dirname(dirname(dirname(d)))), run) ||
        basename(dirname(dirname(d))) != "per_sample_outs" ||
        basename(dirname(dirname(dirname(d)))) != "outs" || basename(d) != "vdj_b")
      fail("[", sid, "] Expected RUN/outs/per_sample_outs/SAMPLE/vdj_b, got ", d)
    if (s$source_project == "GRO11679") {
      expected_run <- sprintf("GRO11679_%s_multi_20260615", s$capture_no)
      if (basename(run) != expected_run || basename(dirname(d)) != s$patient_id)
        fail("[", sid, "] Blood capture/run or hashed sample folder mismatch: ", d)
    } else {
      if (basename(dirname(dirname(run))) != s$capture_no ||
          !grepl("^GRO11735_[a-f0-9]{10}_L5-8$", basename(run)) ||
          basename(dirname(d)) != basename(run))
        fail("[", sid, "] SF capture folder or Cell Ranger run ID mismatch: ", d)
    }
    if (!dir.exists(run) || !dir.exists(d)) fail("[",sid,"] Missing HPC run_dir/vdj_dir: ",d)
    if (normalizePath(d, mustWork=TRUE) !=
        normalizePath(file.path(run, "outs", "per_sample_outs", basename(dirname(d)), "vdj_b"), mustWork=TRUE))
      fail("[", sid, "] Canonical paths inconsistent (possible symlink)")
    # Reject symlinks that redirect into old/foreign alignments.
    if (!startsWith(paste0(normalizePath(d, winslash="/",mustWork=TRUE), "/"), permitted))
      fail("[",sid,"] vdj_dir symlink points outside approved alignment root")
    locations <- file.path(d, needed_files)
    fi <- file.info(locations)
    if (anyNA(fi$size) || any(fi$size <= 0) || any(fi$isdir))
      fail("[",sid,"] Missing/empty VDJ files: ", paste(needed_files[is.na(fi$size) | fi$size <= 0 | (fi$isdir %in% TRUE)], collapse=", "))
    h <- names(read.csv(file.path(d,"filtered_contig_annotations.csv"),nrows=0,check.names=FALSE))
    if (!all(c("barcode","chain") %in% h)) fail("[",sid,"] Unexpected contig CSV schema; inspect header")
    if (!any(c("clonotype_id","raw_clonotype_id") %in% h))
      fail("[",sid,"] No Cell Ranger clonotype ID column")
    if (!any(c("consensus_id","raw_consensus_id") %in% h))
      fail("[",sid,"] No Cell Ranger consensus ID column")
    ch <- names(read.csv(file.path(d,"consensus_annotations.csv"),nrows=0,check.names=FALSE))
    if (!any(c("consensus_id","raw_consensus_id") %in% ch))
      fail("[",sid,"] No consensus ID column in consensus_annotations.csv")
    annotations <- read.csv(file.path(d,"filtered_contig_annotations.csv"),
                            stringsAsFactors=FALSE,check.names=FALSE)
    if (nrow(annotations) == 0L || any(!valid(annotations$barcode)))
      fail("[",sid,"] Empty annotations or missing barcode(s) in Cell Ranger source")
    check_csv_sample_label(annotations,s)
    cap <- s$capture_no[[1L]]
    current_barcodes <- unique(as.character(annotations$barcode))
    previous <- seen_barcodes[[cap]]
    if (!is.null(previous)) for (prior_id in names(previous)) {
      overlap <- intersect(current_barcodes,previous[[prior_id]])
      if (length(overlap)) fail("[",sid,"] ",length(overlap),
          " cell barcodes overlap with ",prior_id," in same capture ",cap,
          ". Refusing potential pooled/duplicated VDJ inputs. Example barcode: ",overlap[[1L]])
    }
    if (is.null(previous)) previous <- list()
    previous[[sid]] <- current_barcodes
    seen_barcodes[[cap]] <- previous
    rm(annotations, current_barcodes)
    say("CHECK OK ",sid," (",s$capture_no,") ",d)
  }
  if (!requireNamespace("Platypus",quietly=TRUE) || !requireNamespace("Biostrings",quietly=TRUE))
    fail("R packages Platypus and Biostrings are required in the active conda environment")
  say("Cohort preflight OK: 16 samples; 13 individuals; 10 captures")
  smp
}

# Compatibility fix for Platypus 3.6.0 with modern Biostrings:
# pairwiseAlignment was moved from Biostrings to pwalign. Avoid recursive
# mutation of the function's language objects; deparse, replace one exact
# *namespace-qualified function name*, and reparse the body. This leaves the
# four requested VDJ_build parameters and the rest of its code unchanged.
make_vdj_build <- function() {
  fn <- get("VDJ_build", asNamespace("Platypus"))
  src <- paste(deparse(body(fn), width.cutoff=500L), collapse="\n")
  needle <- "Biostrings::pairwiseAlignment"
  locations <- gregexpr(needle, src, fixed=TRUE)[[1L]]
  count <- if (locations[[1L]] < 0L) 0L else length(locations)

  if ("pairwiseAlignment" %in% getNamespaceExports("Biostrings")) {
    say("VDJ_build compatibility: unchanged; pairwiseAlignment exported by Biostrings")
    return(list(fn=fn, mode="unmodified", old_calls=count))
  }
  if (count != 1L) fail("Cannot safely adapt Platypus: expected exactly one ",
                        needle, " reference, found ", count)
  if (!requireNamespace("pwalign", quietly=TRUE))
    fail("pwalign package required to replace Biostrings::pairwiseAlignment")

  # Run a real function/API probe before modifying the body.
  probe <- tryCatch(pwalign::pairwiseAlignment(
    pattern=Biostrings::DNAString("ACGTACGT"),
    subject=Biostrings::DNAString("ACGTACGT"),
    type="global-local", gapOpening=10),
    error=function(e) fail("pwalign alignment API probe failed: ", conditionMessage(e)))
  if (!methods::is(probe, "PairwiseAlignments") ||
      !is.character(as.character(probe@subject)))
    fail("pwalign return object differs from the Platypus expectation")

  changed <- sub(needle, "pwalign::pairwiseAlignment", src, fixed=TRUE)
  # Sanity check that undoing our *single replacement* recreates exact source.
  if (!identical(sub("pwalign::pairwiseAlignment", needle, changed, fixed=TRUE), src))
    fail("VDJ_build source changed unexpectedly during compatibility adaptation")
  ast <- tryCatch(parse(text=changed, keep.source=FALSE),
                  error=function(e) fail("Could not parse patched VDJ_build: ",
                                         conditionMessage(e)))
  if (length(ast) != 1L) fail("Unexpected expression count in patched VDJ_build")
  patched <- fn
  body(patched) <- ast[[1L]]
  environment(patched) <- environment(fn)
  say("VDJ_build compatibility: single deparse/parse namespace replacement ",
      "Biostrings -> pwalign; all processing parameters unchanged")
  list(fn=patched, mode="single_source_namespace_replacement_to_pwalign", old_calls=count)
}

# Capture the live R call stack *before* the error unwinds; the previous script
# reported only 'incorrect number of arguments to <-', hiding its origin.
run_vdj_diagnostic <- function(sid, vdj_fn, vdj_dir) {
  withCallingHandlers(
    vdj_fn(
      VDJ.sample.list = c(vdj_dir),
      remove.divergent.cells = TRUE,
      complete.cells.only = TRUE,
      trim.germlines = TRUE
    ),
    error=function(e) {
      dir.create(file.path(results_root, "diagnostics"),
                 showWarnings=FALSE, recursive=TRUE)
      dest <- file.path(results_root, "diagnostics",
                        paste0("ERROR_TRACE_", sid, ".txt"))
      calls <- sys.calls()
      formatted <- vapply(seq_along(calls), function(i) {
        # Deparse with upper bound to avoid dumping long sequence vectors.
        statement <- paste(deparse(calls[[i]], width.cutoff=180L), collapse=" ")
        if (nchar(statement) > 1500L)
          statement <- paste0(substr(statement, 1L, 1500L), " ...[truncated]")
        sprintf("[%03d] %s", i, statement)
      }, character(1))
      lines <- c(paste("Sample:", sid), paste("Time:", Sys.time()),
                 paste("Error:", conditionMessage(e)),
                 paste("Error call:", paste(deparse(conditionCall(e)), collapse=" ")),
                 "\nR call stack (last 75 calls, innermost last):",
                 tail(formatted, 75L))
      try(writeLines(lines, dest), silent=TRUE)
      message("\nVDJ_build failed for ", sid, ": ", conditionMessage(e),
              "\nDiagnostic traceback saved to: ", dest,
              "\nLAST 12 CALLS:\n", paste(tail(formatted, 12L), collapse="\n"))
    }
  )
}

# Cell Ranger v9 headers can use raw_* names where older Platypus expects unprefixed
# identifiers. This temporary, explicitly logged column-alias adapter does not
# alter actual alignments, sequences, barcodes, or ID values.
alias_csv <- function(source, dest) {
  x <- read.csv(source, check.names=FALSE, stringsAsFactors=FALSE)
  for (pair in list(c("raw_clonotype_id","clonotype_id"),c("raw_consensus_id","consensus_id"))) {
    old <- pair[[1L]]; new <- pair[[2L]]
    if (old %in% names(x) && new %in% names(x)) {
      a <- clean_missing(x[[old]]); b <- clean_missing(x[[new]])
      if (!identical(a,b)) fail("Conflicting ",old," and ",new," columns in ",source)
    } else if (old %in% names(x) && !(new %in% names(x))) {
      names(x)[names(x) == old] <- new
      say("ALIAS ",basename(source),": ",old," -> ",new)
    }
  }
  write.csv(x,dest,row.names=FALSE,na="")
  invisible(names(x))
}

prepare_vdj_dir <- function(d, sid) {
  h <- names(read.csv(file.path(d,"filtered_contig_annotations.csv"),nrows=0,check.names=FALSE))
  h2 <- names(read.csv(file.path(d,"consensus_annotations.csv"),nrows=0,check.names=FALSE))
  alias_needed <- any(!c("clonotype_id","consensus_id") %in% h) || !("consensus_id" %in% h2)
  if (!alias_needed) return(list(dir=paste0(d,"/"), adapted=FALSE))
  stage <- file.path(tempdir(),"JIA_CR_schema_stage",sid,"vdj_b")
  if (dir.exists(stage)) fail("Unexpected pre-existing temp schema stage: ",stage)
  dir.create(stage,recursive=TRUE,showWarnings=FALSE)
  for (f in needed_files) {
    source <- file.path(d,f); target <- file.path(stage,f)
    if (f %in% c("filtered_contig_annotations.csv","consensus_annotations.csv")) {
      alias_csv(source,target)
    } else if (!file.symlink(source,target)) {
      if (!file.copy(source,target)) fail("Cannot stage ",source)
    }
  }
  # Keep optional annotation JSON available for older Cell Ranger/Platypus paths.
  extra <- file.path(d,"all_contig_annotations.json")
  if (file.exists(extra)) {
    dest <- file.path(stage,basename(extra))
    if (!file.symlink(extra,dest) && !file.copy(extra,dest)) fail("Cannot stage optional annotation JSON")
  }
  h3 <- names(read.csv(file.path(stage,"filtered_contig_annotations.csv"),nrows=0,check.names=FALSE))
  if (!all(c("barcode","chain","clonotype_id","consensus_id") %in% h3))
    fail("[",sid,"] Cell Ranger schema remains incompatible after aliases")
  say("[",sid,"] USE TEMPORARY SCHEMA ALIASES (original input unchanged): ",stage)
  list(dir=paste0(stage,"/"), adapted=TRUE)
}

atomic_rds <- function(x, path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  if (file.exists(path)) fail("Refusing to overwrite existing output: ",path)
  tmp <- paste0(path,".partial_",Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp),add=TRUE)
  saveRDS(x,tmp,compress="xz")
  if (!file.rename(tmp,path)) fail("Atomic rename failed: ",tmp," -> ",path)
  invisible(path)
}

# Safe union of varying Platypus columns. Retains every column and preserves
# numeric/logical types when compatible; mixed types are converted to character.
combine_two <- function(a,b) {
  if (is.null(a)) return(b)
  columns <- union(names(a),names(b))
  col_type <- function(v) {
    if (is.null(v)) return(NA_character_)
    if (all(is.na(v))) return("missing")
    if (is.factor(v)) return("character")
    ty <- typeof(v)
    if (!ty %in% c("logical","integer","double","character"))
      fail("Unsupported list/complex column type encountered: ",ty)
    ty
  }
  for (nm in columns) {
    av <- if (nm %in% names(a)) a[[nm]] else NULL
    bv <- if (nm %in% names(b)) b[[nm]] else NULL
    ts <- c(col_type(av),col_type(bv))
    ts <- ts[!is.na(ts) & ts != "missing"]
    ty <- if (length(ts)==0L) {
      original <- c(if (!is.null(av)) typeof(av) else NA_character_,
                    if (!is.null(bv)) typeof(bv) else NA_character_)
      original <- original[!is.na(original)]
      if (length(original) && all(original == "logical")) "logical" else "character"
    } else if (all(ts %in% c("integer","double"))) "double"
      else if (length(unique(ts))==1L) ts[[1L]] else "character"
    if (ty == "character" && length(unique(ts)) > 1L)
      say("COLUMN TYPE NOTE: converting heterogeneous column '",nm,"' to character on combine: ",paste(unique(ts),collapse="/"))
    cast <- switch(ty, character=as.character, integer=as.integer,
                   double=as.double, logical=as.logical)
    a[[nm]] <- if (is.null(av)) cast(rep(NA,nrow(a))) else cast(av)
    b[[nm]] <- if (is.null(bv)) cast(rep(NA,nrow(b))) else cast(bv)
  }
  out <- rbind(a[,columns,drop=FALSE],b[,columns,drop=FALSE])
  rownames(out) <- NULL
  out
}

combine_paths <- function(paths) {
  acc <- NULL
  for (p in paths) {
    if (!file.exists(p)) fail("Missing per-sample RDS for combining: ",p)
    one <- readRDS(p)
    acc <- combine_two(acc,one)
    rm(one); invisible(gc(verbose=FALSE))
  }
  acc
}

process_one <- function(s, vdj_fn) {
  sid <- s$sample_id[[1L]]; d <- s$vdj_dir[[1L]]
  say("START sample ",sid," / ",s$capture_no," / ",s$source_project)
  x <- read.csv(file.path(d,"filtered_contig_annotations.csv"),
                stringsAsFactors=FALSE,check.names=FALSE)
  if (nrow(x)==0L || !"barcode" %in% names(x)) fail("Empty/invalid filtered contig CSV")
  input_technical_sample_label <- check_csv_sample_label(x,s)
  if (any(!valid(x$barcode))) fail("Missing source Cell Ranger barcode(s)")
  raw_barcodes <- unique(as.character(x$barcode))
  # Published Platypus removes the literal -1 sequence, not necessarily the
  # original CSV barcode suffix. Verify the mapping is one-to-one.
  key_fn <- function(y) gsub("-1","",as.character(y),fixed=TRUE)
  lookup_keys <- key_fn(raw_barcodes)
  if (anyDuplicated(lookup_keys)) fail("Ambiguous barcode normalization in Cell Ranger input")
  input_clono_col <- if ("raw_clonotype_id" %in% names(x)) "raw_clonotype_id" else "clonotype_id"
  if (!input_clono_col %in% names(x)) fail("Missing source clonotype column")
  raw_clone <- clean_missing(x[[input_clono_col]])
  raw_key <- key_fn(x$barcode)
  clone_sets <- split(raw_clone,raw_key)
  clone_sets <- lapply(clone_sets,function(y) sort(unique(y[valid(y)])))
  prepared <- prepare_vdj_dir(d,sid)
  vdj_dir <- prepared$dir  # trailing slash required by published VDJ_build
  # Same four parameters as v1.2, but retain a useful diagnostic traceback.
  say("[", sid, "] ENTER VDJ_build (remove.divergent.cells=TRUE; ",
      "complete.cells.only=TRUE; trim.germlines=TRUE)")
  z <- run_vdj_diagnostic(sid, vdj_fn, vdj_dir)
  if (!is.data.frame(z) || nrow(z)==0L) fail("VDJ_build returned empty or non-data-frame result")
  if (!"barcode" %in% names(z) || !"clonotype_id" %in% names(z))
    fail("Missing required 'barcode' or 'clonotype_id' in Platypus result")
  if (any(!valid(z$barcode))) fail("Missing barcode from VDJ_build")
  idx <- match(key_fn(z$barcode),lookup_keys)
  if (anyNA(idx)) fail("Platypus barcode(s) not found in Cell Ranger source; cannot preserve raw identity")
  z$barcode_10x_raw <- raw_barcodes[idx]
  z$cell_id <- paste0(sid,"__",z$barcode_10x_raw)
  if (anyDuplicated(z$cell_id)) fail("Multiple VDJ_build rows assigned to one cell; cannot create unique cell_id")
  recovered <- lapply(key_fn(z$barcode),function(k) clone_sets[[k]])
  nclones <- lengths(recovered)
  z$clonotype_id_10x_raw <- vapply(recovered,function(a) if (length(a)==1L) a else NA_character_,"")
  z$clonotype_ids_10x_all_raw <- vapply(recovered,function(a) if (length(a)) paste(a,collapse=";") else NA_character_,"")
  z$clonotype_conflict_10x <- nclones > 1L
  z$clonotype_id_platypus_raw <- clean_missing(z$clonotype_id)
  z$clonotype_id <- NA_character_
  ok <- valid(z$clonotype_id_10x_raw)
  z$clonotype_id[ok] <- paste0(sid,"__",z$clonotype_id_10x_raw[ok])
  z$cellranger_csv_sample_label <- input_technical_sample_label
  z$source_project <- s$source_project; z$capture_no <- s$capture_no
  z$sample_id_platypus_raw <- if ("sample_id" %in% names(z)) as.character(z$sample_id) else NA_character_
  z$sample_id <- sid; z$patient_id <- s$patient_id; z$tissue <- s$tissue; z$side <- s$side
  z$run_dir <- s$run_dir; z$vdj_dir <- d
  z$vdj_dir_build <- vdj_dir; z$schema_alias_applied <- prepared$adapted
  # Validate one row per cell, record paired-chain count only if supported.
  paired <- NA_integer_
  if (all(c("VDJ_chain_count","VJ_chain_count") %in% names(z))) {
    h <- suppressWarnings(as.integer(as.character(z$VDJ_chain_count)))
    l <- suppressWarnings(as.integer(as.character(z$VJ_chain_count)))
    paired <- sum(!is.na(h) & !is.na(l) & h > 0L & l > 0L)
  } else if (all(c("VDJ_chain","VJ_chain") %in% names(z))) {
    paired <- sum(valid(z$VDJ_chain) & valid(z$VJ_chain))
  }
  if (nrow(z) > length(raw_barcodes)) fail("Output unique cells exceed input unique barcode count")
  out <- file.path(results_root,"by_sample",paste0("VDJ_build_",sid,".rds.xz"))
  atomic_rds(z,out)
  ans <- data.frame(source_project=s$source_project,capture_no=s$capture_no,
    sample_id=sid,patient_id=s$patient_id,tissue=s$tissue,side=s$side,
    n_rows=nrow(z),n_unique_cell_barcodes=length(unique(z$barcode_10x_raw)),
    n_unique_clonotypes=length(unique(z$clonotype_id[valid(z$clonotype_id)])),
    n_paired_chain_cells=paired,n_input_barcodes=length(raw_barcodes),
    retention_fraction=nrow(z)/length(raw_barcodes),
    n_ambiguous_input_clonotypes=sum(z$clonotype_conflict_10x),
    cellranger_csv_sample_label=input_technical_sample_label,
    schema_alias_applied=prepared$adapted,run_dir=s$run_dir,vdj_dir=d,
    output_rds=out,stringsAsFactors=FALSE)
  say("SUCCESS sample ",sid,": ",nrow(z)," cells; ",ans$n_unique_clonotypes," clonotypes; saved ",out)
  rm(z,x,raw_clone,recovered,clone_sets); invisible(gc(verbose=FALSE))
  ans
}

smp <- preflight()
starting_mapping_md5 <- unname(tools::md5sum(map_path))
compat <- make_vdj_build()
say("Versions: R=",as.character(getRversion()),"; Platypus=",as.character(packageVersion("Platypus")),
    "; Biostrings=",as.character(packageVersion("Biostrings")),
    "; pwalign=",if(requireNamespace("pwalign",quietly=TRUE)) as.character(packageVersion("pwalign")) else "unavailable")
if (mode == "--check") {
  say("CHECK COMPLETE. No VDJ_build runs or results written.")
  quit(save="no",status=0)
}
if (file.exists(results_root)) fail("RESULTS ALREADY EXIST; refusing to overwrite: ",results_root,
                                    ". For a new run use a new PROJECT_ROOT.")
dir.create(results_root,recursive=TRUE,showWarnings=FALSE)
lock <- file.path(results_root,".RUNNING")
if (!dir.create(lock,showWarnings=FALSE)) fail("Cannot create exclusive run lock: ",lock)
writeLines(c(paste("Started",Sys.time()),paste("Host",Sys.info()[["nodename"]]),
             paste("Mode",mode),paste("Mapping_MD5",unname(tools::md5sum(map_path)))),
           file.path(results_root,"run_provenance.txt"))
if (!file.copy(map_path,file.path(results_root,"sample_mapping_used.tsv"),overwrite=FALSE))
  fail("Could not snapshot mapping TSV")
script_from_args <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if (length(script_from_args)) {
  src <- sub("^--file=","",script_from_args[[1L]])
  file.copy(src,file.path(results_root,"R_source_used.R"),overwrite=FALSE)
}
capture.output({print(sessionInfo()); cat("\nCompatibility mode: ",compat$mode,"\n",sep="")},
               file=file.path(results_root,"sessionInfo.txt"))
write_tsv(smp,file.path(results_root,"sample_input_manifest.tsv"))
status <- tryCatch({
  qc_list <- vector("list",nrow(smp))
  for (i in seq_len(nrow(smp))) {
    sid <- smp$sample_id[[i]]
    qc_list[[i]] <- tryCatch(process_one(smp[i,,drop=FALSE],compat$fn),
      error=function(e) fail("FAILED SAMPLE ",sid,": ",conditionMessage(e)))
  }
  qc <- do.call(rbind,qc_list)
  if (nrow(qc)!=16L || !setequal(qc$sample_id,smp$sample_id)) fail("Internal sample QC mismatch")
  write_tsv(qc,file.path(results_root,"QC_counts_by_sample.tsv"))
  qc_cap <- do.call(rbind,lapply(split(qc,qc$capture_no),function(z) data.frame(
    source_project=unique(z$source_project),capture_no=unique(z$capture_no),
    n_samples=nrow(z),n_patients=length(unique(z$patient_id)),
    n_vdj_build_rows=sum(z$n_rows),
    n_unique_sample_cells=sum(z$n_unique_cell_barcodes),
    n_sample_specific_clonotypes=sum(z$n_unique_clonotypes),
    n_input_sample_barcodes=sum(z$n_input_barcodes),
    n_paired_chain_cells=if(all(is.na(z$n_paired_chain_cells))) NA_integer_ else sum(z$n_paired_chain_cells,na.rm=TRUE),
    stringsAsFactors=FALSE)))
  rownames(qc_cap) <- NULL
  if (nrow(qc_cap)!=10L) fail("Expected 10 capture QC rows")
  write_tsv(qc_cap,file.path(results_root,"QC_counts_by_capture.tsv"))
  for (pid in unique(smp$patient_id)) {
    q <- qc[qc$patient_id==pid,,drop=FALSE]
    p <- combine_paths(q$output_rds)
    if (nrow(p)!=sum(q$n_rows) || any(p$patient_id!=pid) ||
        !setequal(unique(p$sample_id),q$sample_id) || anyDuplicated(p$cell_id))
      fail("Per-patient validation failed: ",pid)
    atomic_rds(p,file.path(results_root,"by_patient",paste0("VDJ_build_",pid,".rds.xz")))
    rm(p); invisible(gc(verbose=FALSE))
    say("SUCCESS patient ",pid,": ",nrow(q)," sample(s)")
  }
  all <- combine_paths(qc$output_rds)
  if (nrow(all)!=sum(qc$n_rows) || length(unique(all$sample_id))!=16L ||
      length(unique(all$patient_id))!=13L || anyDuplicated(all$cell_id) ||
      any(!valid(all$cell_id))) fail("FINAL COMBINED TABLE FAILED VALIDATION")
  if (!identical(unname(tools::md5sum(map_path)), starting_mapping_md5))
    fail("Mapping TSV changed while job was running. Outputs are not marked complete.")
  atomic_rds(all,file.path(results_root,"VDJ_build_all_16_samples.rds.xz"))
  writeLines(c("STATUS: COMPLETE",paste("Finished:",Sys.time()),
    "Input samples: 16","Individuals: 13","Captures: 10",
    paste("VDJ_build rows:",nrow(all)),
    paste("Unique cell IDs:",length(unique(all$cell_id))),
    paste("Sample IDs:",paste(unique(all$sample_id),collapse=", ")),
    paste("Mapping MD5:",unname(tools::md5sum(map_path)))),
    file.path(results_root,"COMPLETION_REPORT.txt"))
  rm(all); invisible(gc(verbose=FALSE))
  TRUE
},error=function(e) {
  msg <- paste("RUN FAILED:",conditionMessage(e))
  message(msg)
  writeLines(c(msg,paste("Time:",Sys.time()),
    "Partial outputs are retained. Do NOT rerun into this results directory.",
    "See diagnostics/ERROR_TRACE_<sample>.txt for the live R call stack."),
    file.path(results_root,"FAILURE_REPORT.txt"))
  FALSE
})
if (!status) quit(save="no",status=1)
if (!file.rename(lock,file.path(results_root,".COMPLETE"))) fail("Could not finalize run lock")
say("PIPELINE COMPLETE: 16/16 samples, 13/13 patients, 10/10 captures")
RSCRIPT

if [[ -e "${R_SCRIPT}" ]]; then
  if ! cmp -s "${TEMP_R}" "${R_SCRIPT}"; then
    echo "ERROR: ${R_SCRIPT} differs from this v1.3 sbatch. Do not mix versions; archive/rename the old generated R script first." >&2
    exit 3
  fi
  rm -f "${TEMP_R}"
else
  mv "${TEMP_R}" "${R_SCRIPT}"
fi
chmod 600 "${R_SCRIPT}"

source "${CONDA_SH}"
conda activate "${R_ENV}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [[ "${MODE}" == "--run" ]]; then
  export TMPDIR="/scratch/${SLURM_JOB_ID:?Must run via sbatch}/jia_vdj_build_v1_3"
else
  export TMPDIR="${PROJECT_ROOT}/01_tmp/preflight_v1_3_${USER:-user}_$$"
fi
mkdir -p "${TMPDIR}"
echo "START: $(date); HOST=$(hostname); MODE=${MODE}; TMPDIR=${TMPDIR}"
Rscript --vanilla "${R_SCRIPT}" "${MODE}" "${PROJECT_ROOT}" "${MAP}"
echo "FINISH: $(date)"
