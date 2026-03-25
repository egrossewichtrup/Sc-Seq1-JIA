#!/usr/bin/env bash
#SBATCH --job-name=AF_build_and_plot_v1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=32:00:00
#SBATCH --output=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AF_pipeline_wo_reclonotyping/%x_%j.out
#SBATCH --error=/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AF_pipeline_wo_reclonotyping/%x_%j.err

set -euo pipefail

###############################################################################
# Paths
###############################################################################
project_root="/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AF_pipeline_wo_reclonotyping"

# Step 1 outputs you already created
input_by_sample="/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/02_results/vdj_build/by_sample"

# Conda and R
conda_sh="/hpc/dla_lti/yermanos_group/software/venvs4/miniforge3/etc/profile.d/conda.sh"
r_env="/hpc/dla_lti/yermanos_group/software/venvs4/r"

###############################################################################
# Shared optional selectors
###############################################################################
ONLY_SAMPLE="${ONLY_SAMPLE:-}"          # example: ONLY_SAMPLE="P001_LK"

###############################################################################
# Step 1 options
###############################################################################
FORCE_REBUILD="${FORCE_REBUILD:-1}"     # set to 1 to overwrite existing outputs
KEEP_SINGLETONS="${KEEP_SINGLETONS:-1}" # set to 1 to keep clonotypes with 1 unique node

###############################################################################
# Step 2 options
###############################################################################
N_CLONES="${N_CLONES:-20}"              # number of trees per sample folder
MIN_NODES="${MIN_NODES:-4}"             # minimum nodes per tree

COLOR_BY="${COLOR_BY:-isotype}"         # isotype, tissue, seurat_clusters
LABEL_BY="${LABEL_BY:-size}"            # name, size, none
EDGE_LABEL="${EDGE_LABEL:-none}"        # original, none
NODE_SIZE="${NODE_SIZE:-expansion}"     # expansion or fixed
SHOW_INNER_NODES="${SHOW_INNER_NODES:-FALSE}"

###############################################################################
# Directories
###############################################################################
mkdir -p "${project_root}/00_config"
mkdir -p "${project_root}/01_logs/slurm"
mkdir -p "${project_root}/02_results/af_objects/by_sample"
mkdir -p "${project_root}/02_results/plots/by_sample"
mkdir -p "${project_root}/02_results/qc"

###############################################################################
# Temp dir
###############################################################################
export TMPDIR="${SLURM_TMPDIR:-/tmp/${USER}/AF_${SLURM_JOB_ID:-manual}}"
mkdir -p "${TMPDIR}"

###############################################################################
# Activate environment
###############################################################################
source "${conda_sh}"
conda activate "${r_env}"

###############################################################################
# Step 1: Build AntibodyForests objects
###############################################################################
export PROJECT_ROOT="${project_root}"
export INPUT_BY_SAMPLE="${input_by_sample}"
export ONLY_SAMPLE
export FORCE_REBUILD
export KEEP_SINGLETONS
export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-1}"

Rscript --vanilla - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(Platypus)
  library(AntibodyForests)
})

options(stringsAsFactors = FALSE)

project_root <- Sys.getenv("PROJECT_ROOT")
input_root   <- Sys.getenv("INPUT_BY_SAMPLE")
only_sample  <- Sys.getenv("ONLY_SAMPLE")
force_rebuild <- as.integer(Sys.getenv("FORCE_REBUILD", "0"))
keep_singletons <- as.integer(Sys.getenv("KEEP_SINGLETONS", "1"))
cpus         <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(cpus) || cpus < 1) cpus <- 1

out_root <- file.path(project_root, "02_results", "af_objects", "by_sample")
qc_root  <- file.path(project_root, "02_results", "qc")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root,  recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_root)) stop("Input by_sample directory not found: ", input_root)

vdj_files <- list.files(
  path = input_root,
  pattern = "^VDJ_build_.*\\.rds\\.xz$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(vdj_files) < 1) stop("No VDJ_build_*.rds.xz files found under: ", input_root)

get_sample_id <- function(p) {
  b <- basename(p)
  b <- sub("^VDJ_build_", "", b)
  b <- sub("\\.rds\\.xz$", "", b)
  b
}

sample_ids <- vapply(vdj_files, get_sample_id, character(1))
o <- order(sample_ids)
vdj_files <- vdj_files[o]
sample_ids <- sample_ids[o]

if (nzchar(only_sample)) {
  keep <- sample_ids == only_sample
  if (!any(keep)) stop("ONLY_SAMPLE was set but not found: ", only_sample)
  vdj_files <- vdj_files[keep]
  sample_ids <- sample_ids[keep]
}

manifest_path <- file.path(project_root, "00_config", "samples_for_af_build.tsv")
write.table(
  data.frame(sample_id = sample_ids, vdj_rds = vdj_files),
  manifest_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

qc_path <- file.path(qc_root, "Af_build_summary_by_sample.tsv")
qc_header <- data.frame(
  sample_id = character(0),
  input_rds = character(0),
  status = character(0),
  mode = character(0),
  clonotype_source = character(0),
  n_rows_input = integer(0),
  n_rows_after_rowfilter = integer(0),
  n_clonotypes_total = integer(0),
  n_clonotypes_kept = integer(0),
  n_rows_final = integer(0),
  seconds = numeric(0),
  error = character(0),
  stringsAsFactors = FALSE
)
write.table(qc_header, qc_path, sep = "\t", quote = FALSE, row.names = FALSE)

append_qc <- function(row_df) {
  write.table(row_df, qc_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

safe_nzchar <- function(x) {
  !is.na(x) && nchar(as.character(x)) > 0
}

pick_node_features <- function(vdj) {
  candidates <- c("isotype", "tissue", "side", "seurat_clusters", "sample_id", "patient_id", "capture_id")
  intersect(candidates, colnames(vdj))
}

choose_mode <- function(vdj) {
  modes <- list(
    list(name="heavy_plus_light_aa_trimmed",
         seq=c("VDJ_sequence_aa_trimmed","VJ_sequence_aa_trimmed"),
         gl =c("VDJ_germline_aa_trimmed","VJ_germline_aa_trimmed")),
    list(name="heavy_plus_light_nt_trimmed",
         seq=c("VDJ_sequence_nt_trimmed","VJ_sequence_nt_trimmed"),
         gl =c("VDJ_germline_nt_trimmed","VJ_germline_nt_trimmed")),
    list(name="heavy_only_aa_trimmed",
         seq=c("VDJ_sequence_aa_trimmed"),
         gl =c("VDJ_germline_aa_trimmed")),
    list(name="heavy_only_nt_trimmed",
         seq=c("VDJ_sequence_nt_trimmed"),
         gl =c("VDJ_germline_nt_trimmed"))
  )

  for (m in modes) {
    if (all(m$seq %in% colnames(vdj)) && all(m$gl %in% colnames(vdj))) {
      ok_rows <- rep(TRUE, nrow(vdj))
      for (k in seq_along(m$seq)) {
        ok_rows <- ok_rows & vapply(vdj[[m$seq[k]]], safe_nzchar, logical(1)) &
          vapply(vdj[[m$gl[k]]], safe_nzchar, logical(1))
      }
      if (sum(ok_rows) >= 1) return(m)
    }
  }
  NULL
}

pick_safe_clonotype <- function(vdj) {
  if ("clonotype_id_10x_raw" %in% colnames(vdj)) return("clonotype_id_10x_raw")
  if ("clonotype_id_10x" %in% colnames(vdj)) return("clonotype_id_10x")
  if ("clonotype_id" %in% colnames(vdj)) return("clonotype_id")
  if ("clonotype" %in% colnames(vdj)) return("clonotype")
  NA_character_
}

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  in_rds <- vdj_files[i]

  sample_dir <- file.path(out_root, sid)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  out_rds     <- file.path(sample_dir, paste0("AntibodyForests_", sid, ".rds.xz"))
  out_session <- file.path(sample_dir, paste0("AntibodyForests_", sid, ".sessionInfo.txt"))
  out_log     <- file.path(sample_dir, paste0("AntibodyForests_", sid, ".run.log"))
  out_err     <- file.path(sample_dir, paste0("AntibodyForests_", sid, ".error.txt"))

  if (file.exists(out_rds) && force_rebuild == 0) {
    info <- file.info(out_rds)
    if (!is.na(info$size) && info$size > 0) {
      append_qc(data.frame(
        sample_id = sid,
        input_rds = in_rds,
        status = "skipped_existing_output",
        mode = NA_character_,
        clonotype_source = NA_character_,
        n_rows_input = NA_integer_,
        n_rows_after_rowfilter = NA_integer_,
        n_clonotypes_total = NA_integer_,
        n_clonotypes_kept = NA_integer_,
        n_rows_final = NA_integer_,
        seconds = 0,
        error = NA_character_,
        stringsAsFactors = FALSE
      ))
      next
    }
  }

  t0 <- Sys.time()

  res_row <- data.frame(
    sample_id = sid,
    input_rds = in_rds,
    status = "error",
    mode = NA_character_,
    clonotype_source = NA_character_,
    n_rows_input = NA_integer_,
    n_rows_after_rowfilter = NA_integer_,
    n_clonotypes_total = NA_integer_,
    n_clonotypes_kept = NA_integer_,
    n_rows_final = NA_integer_,
    seconds = NA_real_,
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  con <- file(out_log, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")

  tryCatch({

    vdj <- readRDS(in_rds)
    res_row$n_rows_input <- nrow(vdj)

    if (!("barcode" %in% colnames(vdj))) stop("Missing required column: barcode")

    if (!("sample_id" %in% colnames(vdj))) {
      vdj$sample_id <- sid
      message("Inserted sample_id column with value: ", sid)
    }

    cl_source <- pick_safe_clonotype(vdj)
    if (is.na(cl_source)) stop("No clonotype column found. Expected clonotype_id_10x_raw or clonotype_id.")
    res_row$clonotype_source <- cl_source

    vdj$clonotype_id <- as.character(vdj[[cl_source]])

    if (any(grepl("_", vdj$clonotype_id, fixed = TRUE), na.rm = TRUE)) {
      vdj$clonotype_id_original <- vdj$clonotype_id
      vdj$clonotype_id <- gsub("_", ".", vdj$clonotype_id, fixed = TRUE)
      message("Sanitized underscores in clonotype_id by replacing with '.'")
    }

    vdj <- vdj[vapply(vdj$clonotype_id, safe_nzchar, logical(1)), , drop = FALSE]
    if (nrow(vdj) < 1) {
      res_row$status <- "skipped_no_clonotypes"
      res_row$error <- "All clonotype ids were empty or NA after filtering."
      res_row$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      append_qc(res_row)
      next
    }

    mode <- choose_mode(vdj)
    if (is.null(mode)) stop("No usable sequence and germline columns found in this VDJ object.")
    res_row$mode <- mode$name

    keep_rows <- rep(TRUE, nrow(vdj))
    for (k in seq_along(mode$seq)) {
      keep_rows <- keep_rows &
        vapply(vdj[[mode$seq[k]]], safe_nzchar, logical(1)) &
        vapply(vdj[[mode$gl[k]]], safe_nzchar, logical(1))
    }
    vdj2 <- vdj[keep_rows, , drop = FALSE]
    res_row$n_rows_after_rowfilter <- nrow(vdj2)

    if (nrow(vdj2) < 2) {
      res_row$status <- "skipped_too_few_rows_after_filter"
      res_row$error <- "Fewer than 2 rows remain after sequence and germline filtering."
      res_row$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      append_qc(res_row)
      next
    }

    node_key <- apply(vdj2[, mode$seq, drop = FALSE], 1, function(z) paste(z, collapse="|"))
    clono_ids <- vdj2$clonotype_id
    clono_nodes <- tapply(node_key, clono_ids, function(x) length(unique(x)))

    res_row$n_clonotypes_total <- length(clono_nodes)

    if (keep_singletons == 1) {
      keep_clonos <- names(clono_nodes)
    } else {
      keep_clonos <- names(clono_nodes)[clono_nodes >= 2]
    }
    res_row$n_clonotypes_kept <- length(keep_clonos)

    vdj3 <- vdj2[vdj2$clonotype_id %in% keep_clonos, , drop = FALSE]
    res_row$n_rows_final <- nrow(vdj3)

    if (nrow(vdj3) < 2 || length(keep_clonos) < 1) {
      res_row$status <- "skipped_no_clonotypes_kept"
      res_row$error <- "No clonotypes remained after filtering."
      res_row$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      append_qc(res_row)
      next
    }

    node_features <- pick_node_features(vdj3)

    af <- AntibodyForests::Af_build(
      VDJ = vdj3,
      construction.method = "phylo.network.default",
      sequence.columns = if (length(mode$seq) == 1) mode$seq[[1]] else mode$seq,
      germline.columns = if (length(mode$gl) == 1) mode$gl[[1]] else mode$gl,
      concatenate.sequences = FALSE,
      node.features = node_features,
      parallel = FALSE
    )

    saveRDS(af, out_rds, compress = "xz")
    writeLines(capture.output(sessionInfo()), out_session)

    res_row$status <- "ok"
    res_row$error <- NA_character_
    res_row$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    append_qc(res_row)

    rm(vdj, vdj2, vdj3, af)
    gc()

  }, error = function(e) {
    msg <- conditionMessage(e)
    res_row$status <- "error"
    res_row$error <- msg
    res_row$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    calls <- sys.calls()
    call_lines <- vapply(seq_along(calls), function(ii) {
      paste0(ii, ": ", paste(deparse(calls[[ii]]), collapse = " "))
    }, character(1))

    writeLines(
      c(
        paste0("Sample: ", sid),
        paste0("Input: ", in_rds),
        paste0("Error: ", msg),
        "",
        "Call stack:",
        call_lines
      ),
      out_err
    )

    append_qc(res_row)
  })

  sink(type = "message")
  sink(type = "output")
  close(con)
}

cat("Wrote manifest: ", manifest_path, "\n", sep = "")
cat("Wrote QC summary: ", qc_path, "\n", sep = "")
cat("Done\n")
RSCRIPT

###############################################################################
# Step 2: Plot trees from the built AntibodyForests objects
###############################################################################
export PROJECT_ROOT="${project_root}"
export ONLY_SAMPLE
export N_CLONES
export MIN_NODES
export COLOR_BY
export LABEL_BY
export EDGE_LABEL
export NODE_SIZE
export SHOW_INNER_NODES

Rscript --vanilla - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(AntibodyForests)
  library(igraph)
})

options(stringsAsFactors = FALSE)

project_root <- Sys.getenv("PROJECT_ROOT")
only_sample  <- Sys.getenv("ONLY_SAMPLE")

n_clones  <- as.integer(Sys.getenv("N_CLONES", "20"))
min_nodes <- as.integer(Sys.getenv("MIN_NODES", "4"))
if (is.na(n_clones) || n_clones < 1) n_clones <- 20
if (is.na(min_nodes) || min_nodes < 1) min_nodes <- 4

color_by <- Sys.getenv("COLOR_BY", "isotype")
label_by <- Sys.getenv("LABEL_BY", "size")
edge_label <- Sys.getenv("EDGE_LABEL", "none")
node_size <- Sys.getenv("NODE_SIZE", "expansion")
show_inner_nodes <- tolower(Sys.getenv("SHOW_INNER_NODES", "FALSE")) %in% c("true","t","1","yes","y")

# Second PDF node labels: node identifier (igraph vertex name)
# If you want sample labels instead, change "name" to "sample_id" here.
label_by_nodeid <- "name"

af_root   <- file.path(project_root, "02_results", "af_objects", "by_sample")
plot_root <- file.path(project_root, "02_results", "plots", "by_sample")
qc_root   <- file.path(project_root, "02_results", "qc")
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root,   recursive = TRUE, showWarnings = FALSE)

af_files <- list.files(
  path = af_root,
  pattern = "^AntibodyForests_.*\\.rds\\.xz$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(af_files) < 1) stop("No AntibodyForests_*.rds.xz files found under: ", af_root)

sid_from_path <- function(p) basename(dirname(p))
sids <- vapply(af_files, sid_from_path, character(1))
o <- order(sids)
af_files <- af_files[o]
sids <- sids[o]

if (nzchar(only_sample)) {
  keep <- sids == only_sample
  if (!any(keep)) stop("ONLY_SAMPLE was set, but not found: ", only_sample)
  af_files <- af_files[keep]
  sids <- sids[keep]
}

split_key <- function(key) {
  key <- as.character(key)
  m1 <- regexec("^(.+)\\.(clonotype[^.]+)$", key)
  r1 <- regmatches(key, m1)[[1]]
  if (length(r1) == 3) return(list(sample = r1[2], clonotype = r1[3], key = key))

  m2 <- regexec("^(.+)__(clonotype[^.]+)$", key)
  r2 <- regmatches(key, m2)[[1]]
  if (length(r2) == 3) return(list(sample = r2[2], clonotype = r2[3], key = key))

  list(sample = NA_character_, clonotype = NA_character_, key = key)
}

collect_graphs_recursive <- function(x, cur_name = NA_character_, out = NULL) {
  if (is.null(out)) out <- list()

  if (inherits(x, "igraph")) {
    out[[length(out) + 1]] <- list(key = cur_name, graph = x)
    return(out)
  }

  if (is.list(x)) {
    nms <- names(x)

    if (!is.null(nms) && "igraph" %in% nms && inherits(x[["igraph"]], "igraph")) {
      out[[length(out) + 1]] <- list(key = cur_name, graph = x[["igraph"]])
    }

    if (!is.null(nms)) {
      for (nm in nms) out <- collect_graphs_recursive(x[[nm]], cur_name = nm, out = out)
    } else {
      for (k in seq_along(x)) out <- collect_graphs_recursive(x[[k]], cur_name = as.character(k), out = out)
    }
  }

  out
}

call_Af_plot_tree <- function(obj, sample_name, clonotype_name,
                              color_by, label_by, edge_label,
                              node_size, show_inner_nodes,
                              main_title, sub_title) {
  fml <- names(formals(AntibodyForests::Af_plot_tree))
  args <- list(
    AntibodyForests_object = obj,
    sample = sample_name,
    clonotype = clonotype_name
  )
  if ("color.by" %in% fml) args$color.by <- color_by
  if ("label.by" %in% fml) args$label.by <- label_by
  if ("edge.label" %in% fml) args$edge.label <- edge_label
  if ("node.size" %in% fml) args$node.size <- node_size
  if ("show.inner.nodes" %in% fml) args$show.inner.nodes <- show_inner_nodes
  if ("main.title" %in% fml) args$main.title <- main_title
  if ("sub.title" %in% fml) args$sub.title <- sub_title
  do.call(AntibodyForests::Af_plot_tree, args)
}

plot_igraph_fallback <- function(g, color_by, label_by, edge_label, node_size, main_title, sub_title) {
  vattrs <- vertex_attr_names(g)

  lab <- rep("", vcount(g))
  if (label_by == "name" && "name" %in% vattrs) lab <- as.character(V(g)$name)
  if (label_by == "size" && "size" %in% vattrs) lab <- as.character(V(g)$size)
  if (label_by == "size" && !("size" %in% vattrs) && "expansion" %in% vattrs) lab <- as.character(V(g)$expansion)
  if (label_by == "sample_id" && "sample_id" %in% vattrs) lab <- as.character(V(g)$sample_id)

  vsz <- rep(12, vcount(g))
  if (node_size == "expansion") {
    if ("size" %in% vattrs) vsz <- pmax(6, pmin(40, as.numeric(V(g)$size)))
    if (!("size" %in% vattrs) && "expansion" %in% vattrs) vsz <- pmax(6, pmin(40, as.numeric(V(g)$expansion)))
  }

  vcol <- rep("grey70", vcount(g))
  if (color_by %in% vattrs) {
    vals <- as.character(vertex_attr(g, color_by))
    lev <- sort(unique(vals))
    pal <- grDevices::rainbow(length(lev))
    vcol <- pal[match(vals, lev)]
    legend("topleft", legend = lev, col = pal, pch = 19, cex = 0.7, bty = "n")
  }

  elab <- NA
  if (edge_label == "original") {
    if ("weight" %in% edge_attr_names(g)) elab <- as.character(E(g)$weight)
    if (!("weight" %in% edge_attr_names(g)) && "distance" %in% edge_attr_names(g)) elab <- as.character(E(g)$distance)
  }

  plot(
    g,
    vertex.label = lab,
    vertex.size = vsz,
    vertex.color = vcol,
    edge.label = if (all(is.na(elab))) NA else elab,
    main = main_title
  )
  mtext(sub_title, side = 3, line = 0.2, cex = 0.8)
}

qc_rows <- list()

for (i in seq_along(sids)) {
  sid <- sids[[i]]
  in_rds <- af_files[[i]]

  out_dir <- file.path(plot_root, sid)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  run_log <- file.path(out_dir, paste0("Af_plot_tree_", sid, ".run.log"))
  rank_tsv <- file.path(out_dir, paste0("clonotype_ranking_", sid, ".tsv"))
  pdf_out <- file.path(out_dir, paste0("trees_top", n_clones, "_minNodes", min_nodes, "_color", color_by, "_label", label_by, ".pdf"))

  # Second PDF with node identifiers
  pdf_out_nodeid <- file.path(out_dir, paste0("trees_top", n_clones, "_minNodes", min_nodes, "_color", color_by, "_labelnodeid.pdf"))

  if (file.exists(pdf_out)) file.remove(pdf_out)
  if (file.exists(pdf_out_nodeid)) file.remove(pdf_out_nodeid)

  t0 <- Sys.time()

  log_con <- file(run_log, open = "wt")
  sink(log_con, type = "output")
  sink(log_con, type = "message")

  cat("Sample folder: ", sid, "\n", sep = "")
  cat("Input AF RDS: ", in_rds, "\n", sep = "")
  cat("Output PDF:  ", pdf_out, "\n", sep = "")
  cat("Output PDF (node ids):  ", pdf_out_nodeid, "\n", sep = "")
  cat("Settings: n_clones=", n_clones,
      " min_nodes=", min_nodes,
      " color_by=", color_by,
      " label_by=", label_by,
      " edge_label=", edge_label,
      " node_size=", node_size,
      " show_inner_nodes=", show_inner_nodes,
      "\n", sep = "")

  status <- "error"
  err_msg <- NA_character_
  n_found <- 0L
  n_ranked <- 0L
  n_plotted <- 0L
  n_plotted_af <- 0L
  n_plotted_fallback <- 0L

  tryCatch({
    af <- readRDS(in_rds)
    af_list <- list(Default = af)

    rec <- collect_graphs_recursive(af, cur_name = NA_character_)
    if (length(rec) < 1) stop("No igraph objects found anywhere inside this AntibodyForests RDS.")

    keys <- vapply(rec, function(z) z$key, character(1))
    ok_key <- !is.na(keys) & nzchar(keys) & keys != "igraph"
    rec <- rec[ok_key]
    if (length(rec) < 1) stop("Igraph objects were found, but no usable keys were associated with them.")

    keys <- vapply(rec, function(z) z$key, character(1))
    rec <- rec[!duplicated(keys)]
    keys <- vapply(rec, function(z) z$key, character(1))
    n_found <- length(rec)

    info <- lapply(seq_along(rec), function(j) {
      g <- rec[[j]]$graph
      k <- rec[[j]]$key
      sp <- split_key(k)
      data.frame(
        key = k,
        parsed_sample = sp$sample,
        parsed_clonotype = sp$clonotype,
        n_nodes = igraph::vcount(g),
        n_edges = igraph::ecount(g),
        stringsAsFactors = FALSE
      )
    })
    info <- do.call(rbind, info)
    info <- info[!is.na(info$parsed_sample) & !is.na(info$parsed_clonotype), , drop = FALSE]
    if (nrow(info) < 1) stop("Trees were found, but none had parseable keys with a clonotype suffix.")

    info <- info[info$n_nodes >= min_nodes, , drop = FALSE]
    info <- info[order(-info$n_nodes, -info$n_edges), , drop = FALSE]
    write.table(info, rank_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    n_ranked <- nrow(info)

    if (nrow(info) < 1) {
      status <- "skipped_no_trees_passing_min_nodes"
      err_msg <- paste0("No trees with at least ", min_nodes, " nodes.")
    } else {
      top <- info[seq_len(min(nrow(info), n_clones)), , drop = FALSE]

      # PDF 1: unchanged behavior
      pdf(pdf_out, width = 10, height = 8, onefile = TRUE, useDingbats = FALSE)
      on.exit(dev.off(), add = TRUE)

      for (j in seq_len(nrow(top))) {
        k <- top$key[[j]]
        s_parsed <- top$parsed_sample[[j]]
        c_parsed <- top$parsed_clonotype[[j]]

        g <- rec[[match(k, keys)]]$graph

        main_title <- paste0("Lineage tree, folder ", sid)
        sub_title <- paste0("key=", k, " nodes=", top$n_nodes[[j]])

        cat("Plotting: key=", k, " parsed_sample=", s_parsed, " parsed_clonotype=", c_parsed, "\n", sep = "")

        sample_candidates <- unique(c(s_parsed, sid))
        clonotype_candidates <- unique(c(c_parsed, k, paste0(s_parsed, ".", c_parsed), paste0(s_parsed, "__", c_parsed)))

        plotted_this <- FALSE
        last_err <- NA_character_

        for (s_try in sample_candidates) {
          for (c_try in clonotype_candidates) {
            for (obj_try in list(af_list, af)) {
              tryCatch({
                call_Af_plot_tree(
                  obj = obj_try,
                  sample_name = s_try,
                  clonotype_name = c_try,
                  color_by = color_by,
                  label_by = label_by,
                  edge_label = edge_label,
                  node_size = node_size,
                  show_inner_nodes = show_inner_nodes,
                  main_title = main_title,
                  sub_title = paste0("Af_plot_tree: sample=", s_try, " clonotype=", c_try)
                )
                plotted_this <- TRUE
                n_plotted <- n_plotted + 1L
                n_plotted_af <- n_plotted_af + 1L
                cat("SUCCESS Af_plot_tree with sample=", s_try, " clonotype=", c_try, "\n", sep = "")
                break
              }, error = function(e2) {
                last_err <<- conditionMessage(e2)
              })
              if (plotted_this) break
            }
            if (plotted_this) break
          }
          if (plotted_this) break
        }

        if (!plotted_this) {
          cat("Af_plot_tree failed for key ", k, ". Last error: ", last_err, "\n", sep = "")
          cat("Falling back to direct igraph plotting for key ", k, "\n", sep = "")
          plot_igraph_fallback(
            g = g,
            color_by = color_by,
            label_by = label_by,
            edge_label = edge_label,
            node_size = node_size,
            main_title = main_title,
            sub_title = paste0("Fallback igraph plot: ", sub_title)
          )
          n_plotted <- n_plotted + 1L
          n_plotted_fallback <- n_plotted_fallback + 1L
        }
      }

      # PDF 2: node identifiers
      pdf(pdf_out_nodeid, width = 10, height = 8, onefile = TRUE, useDingbats = FALSE)
      on.exit(dev.off(), add = TRUE)

      for (j in seq_len(nrow(top))) {
        k <- top$key[[j]]
        s_parsed <- top$parsed_sample[[j]]
        c_parsed <- top$parsed_clonotype[[j]]

        g <- rec[[match(k, keys)]]$graph

        main_title <- paste0("Lineage tree, folder ", sid)
        sub_title <- paste0("key=", k, " nodes=", top$n_nodes[[j]])

        cat("Plotting node-id PDF: key=", k, " parsed_sample=", s_parsed, " parsed_clonotype=", c_parsed, "\n", sep = "")

        sample_candidates <- unique(c(s_parsed, sid))
        clonotype_candidates <- unique(c(c_parsed, k, paste0(s_parsed, ".", c_parsed), paste0(s_parsed, "__", c_parsed)))

        plotted_this <- FALSE
        last_err <- NA_character_

        for (s_try in sample_candidates) {
          for (c_try in clonotype_candidates) {
            for (obj_try in list(af_list, af)) {
              tryCatch({
                call_Af_plot_tree(
                  obj = obj_try,
                  sample_name = s_try,
                  clonotype_name = c_try,
                  color_by = color_by,
                  label_by = label_by_nodeid,
                  edge_label = edge_label,
                  node_size = node_size,
                  show_inner_nodes = show_inner_nodes,
                  main_title = main_title,
                  sub_title = paste0("Af_plot_tree node-id: sample=", s_try, " clonotype=", c_try)
                )
                plotted_this <- TRUE
                cat("SUCCESS Af_plot_tree node-id with sample=", s_try, " clonotype=", c_try, "\n", sep = "")
                break
              }, error = function(e2) {
                last_err <<- conditionMessage(e2)
              })
              if (plotted_this) break
            }
            if (plotted_this) break
          }
          if (plotted_this) break
        }

        if (!plotted_this) {
          cat("Af_plot_tree node-id failed for key ", k, ". Last error: ", last_err, "\n", sep = "")
          cat("Falling back to direct igraph plotting (node-id) for key ", k, "\n", sep = "")
          plot_igraph_fallback(
            g = g,
            color_by = color_by,
            label_by = label_by_nodeid,
            edge_label = edge_label,
            node_size = node_size,
            main_title = main_title,
            sub_title = paste0("Fallback igraph plot node-id: ", sub_title)
          )
        }
      }

      if (n_plotted > 0) {
        status <- "ok"
        err_msg <- NA_character_
      } else {
        status <- "error_no_plots_produced"
        err_msg <- "No plots were produced."
      }
    }

    rm(af)
    gc()

  }, error = function(e) {
    err_msg <<- conditionMessage(e)
    status <<- "error"
    cat("ERROR: ", err_msg, "\n", sep = "")
  })

  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  pdf_size <- if (file.exists(pdf_out)) file.info(pdf_out)$size else 0
  pdf_size_nodeid <- if (file.exists(pdf_out_nodeid)) file.info(pdf_out_nodeid)$size else 0

  sink(type = "message")
  sink(type = "output")
  close(log_con)

  qc_rows[[length(qc_rows) + 1]] <- data.frame(
    sample_id = sid,
    input_rds = in_rds,
    status = status,
    n_graphs_found = n_found,
    n_trees_ranked_ge_min_nodes = n_ranked,
    n_plotted_total = n_plotted,
    n_plotted_Af_plot_tree = n_plotted_af,
    n_plotted_fallback_igraph = n_plotted_fallback,
    n_clones_requested = n_clones,
    min_nodes = min_nodes,
    color_by = color_by,
    label_by = label_by,
    edge_label = edge_label,
    node_size = node_size,
    show_inner_nodes = show_inner_nodes,
    pdf_bytes = as.numeric(pdf_size),
    pdf_bytes_nodeid = as.numeric(pdf_size_nodeid),
    seconds = secs,
    error = err_msg,
    stringsAsFactors = FALSE
  )
}

qc_df <- do.call(rbind, qc_rows)
qc_path <- file.path(qc_root, "Af_plot_tree_summary_by_sample.tsv")
write.table(qc_df, qc_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote plot QC summary: ", qc_path, "\n", sep = "")
RSCRIPT

###############################################################################
# Step 3: Merge per-sample PDFs into 2 combined PDFs (no R dependencies)
###############################################################################
plot_root="${project_root}/02_results/plots/by_sample"

# Try to make a merge tool available via modules (harmless if "module" not present)
if command -v module >/dev/null 2>&1; then
  module --silent avail qpdf ghostscript poppler >/dev/null 2>&1 || true
  module --silent load qpdf 2>/dev/null || true
  module --silent load ghostscript 2>/dev/null || true
  module --silent load poppler 2>/dev/null || true   # provides pdfunite sometimes
fi

merge_one_type () {
  local pattern="$1"
  local outpdf="$2"

  cd "$plot_root"

  mapfile -d '' files < <(find . -type f -name "$pattern" -print0 | sort -z)

  if [[ "${#files[@]}" -eq 0 ]]; then
    echo "No files found for pattern: $pattern"
    return 0
  fi

  echo "Merging ${#files[@]} PDFs into: ${plot_root}/${outpdf}"

  if command -v pdfunite >/dev/null 2>&1; then
    pdfunite "${files[@]}" "${outpdf}"
    return $?
  fi

  if command -v qpdf >/dev/null 2>&1; then
    qpdf --empty --pages "${files[@]}" -- "${outpdf}"
    return $?
  fi

  if command -v gs >/dev/null 2>&1; then
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="${outpdf}" "${files[@]}"
    return $?
  fi

  echo "ERROR: No PDF merge tool found (pdfunite, qpdf, gs)."
  echo "Hint: check modules: module avail qpdf ghostscript poppler"
  return 1
}

# Merge node-id PDFs
MERGE_PATTERN="trees_top${N_CLONES}_minNodes${MIN_NODES}_color${COLOR_BY}_labelnodeid.pdf"
merge_one_type "${MERGE_PATTERN}" "merged_${MERGE_PATTERN}"

# Merge size PDFs (your main output)
MERGE_PATTERN="trees_top${N_CLONES}_minNodes${MIN_NODES}_color${COLOR_BY}_label${LABEL_BY}.pdf"
merge_one_type "${MERGE_PATTERN}" "merged_${MERGE_PATTERN}"

echo "Finished"

