
# create_combined_barplot_HS_vs_genetic_demultiplexing — README
#
# • What this sbatch script does
#   - Builds a per-capture, per-sample cell-count comparison between:
#       1) HashSolo demultiplexing (from the demux pipeline’s STEP 6 output table)
#       2) Genetic demultiplexing (Souporcell) after “patient/sample” annotation
#   - Writes two tables + a combined (dodged) barplot figure (PNG + PDF).
#
# • Inputs (auto-discovered)
#   - HashSolo counts table (STEP 6):
#     /hpc/dla_lti/yermanos_group/GRO11679/demultiplexing/seurat/<TS>_step4_3_split/indices/counts_per_capture_sample.csv
#     Columns expected: HS___capture__, HS_Sample_short, n_cells
#
#   - Genetic demux annotated tables:
#     /hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing/demux/**/souporcell_with_patient_*/<CAPTURE>/<CAPTURE>_clusters_with_patient.tsv
#     Columns used: status, patient_id_from_HTO
#
# • How counts are computed
#   - HashSolo:
#       - Read counts_per_capture_sample.csv as-is and standardize to: capture, sample_short, n_cells
#   - Genetic demultiplexing:
#       - For each *_clusters_with_patient.tsv:
#           - keep rows where status == "singlet" and patient_id_from_HTO is non-empty
#           - count rows per (capture, patient_id_from_HTO)
#
# • How the plot is made
#   - Combine both methods into one long table (capture, sample_short, method, n_cells)
#   - Expand methods only *within* existing (capture, sample_short) pairs so each capture facet
#     shows only its own samples (no cross-capture zero-filled samples).
#   - Plot as dodged bars (HashSolo vs Genetic) and facet by capture.
#
# • Outputs (written to a timestamped OUTDIR)
#   - hashsolo_counts_per_capture_sample.tsv
#   - genetic_counts_per_capture_sample.tsv
#   - combined_counts_long.tsv
#   - HS_vs_genetic_counts_per_capture.png
#   - HS_vs_genetic_counts_per_capture.pdf












#!/usr/bin/env bash
#SBATCH --job-name=HS_vs_GenDem_barplot
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --chdir=/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing
#SBATCH --output=/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing/logs/%x_%j.out
#SBATCH --error=/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing/logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=egrossewichtrup@umcutrecht.nl

set -euo pipefail
log(){ echo "[$(date +%F_%T)] $*"; }

PROJ_MAIN="/hpc/dla_lti/yermanos_group/GRO11679"
DEMUX_ROOT="${PROJ_MAIN}/demultiplexing"
GENETIC_ROOT="${PROJ_MAIN}/genetic_demultiplexing"

R_BIN="/hpc/dla_lti/yermanos_group/software/venvs4/r/bin/Rscript"
if [[ ! -x "$R_BIN" ]]; then
  R_BIN="$(command -v Rscript)"
fi

mkdir -p "${GENETIC_ROOT}/logs" "${GENETIC_ROOT}/tmp"

# 1) Find HashSolo counts table (STEP 6 split indices output)
HS_COUNTS_CSV="$(find "${DEMUX_ROOT}/seurat" -type f -path "*_step4_3_split/indices/counts_per_capture_sample.csv" 2>/dev/null | sort | tail -n1 || true)"
if [[ -z "${HS_COUNTS_CSV}" || ! -s "${HS_COUNTS_CSV}" ]]; then
  log "FATAL: Could not find HashSolo counts_per_capture_sample.csv under ${DEMUX_ROOT}/seurat/*_step4_3_split/indices/"
  exit 2
fi

# 2) Confirm genetic demux annotated tables exist
GEN_FILES_N="$(find "${GENETIC_ROOT}/demux" -type f -name "*_clusters_with_patient.tsv" 2>/dev/null | wc -l | tr -d ' ')"
if [[ "${GEN_FILES_N}" -le 0 ]]; then
  log "FATAL: No genetic demux *_clusters_with_patient.tsv found under ${GENETIC_ROOT}/demux"
  exit 3
fi

TS="$(date +%Y%m%d_%H%M%S)"
OUTDIR="${GENETIC_ROOT}/combined_HS_vs_genetic_${TS}"
mkdir -p "$OUTDIR"

JOBTMP="${GENETIC_ROOT}/tmp/HS_vs_GenDem_${SLURM_JOB_ID:-noid}_${TS}"
mkdir -p "$JOBTMP"
export TMPDIR="$JOBTMP" TMP="$JOBTMP" TEMP="$JOBTMP"

log "Using:"
log "  HS_COUNTS_CSV = $HS_COUNTS_CSV"
log "  GENETIC_ROOT  = $GENETIC_ROOT"
log "  OUTDIR        = $OUTDIR"
log "  TMPDIR        = $TMPDIR"
log "  R_BIN         = $R_BIN"

export HS_COUNTS_CSV GENETIC_ROOT OUTDIR

"$R_BIN" - <<'RS'
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), paste0(...)))

hs_path <- Sys.getenv("HS_COUNTS_CSV")
gen_root <- Sys.getenv("GENETIC_ROOT")
outdir <- Sys.getenv("OUTDIR")

stopifnot(file.exists(hs_path), dir.exists(gen_root), outdir != "")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 1) Read HashSolo counts table (STEP 6 indices)
# -------------------------
hs <- readr::read_csv(hs_path, show_col_types = FALSE)
req_hs <- c("HS___capture__", "HS_Sample_short", "n_cells")
miss_hs <- setdiff(req_hs, names(hs))
if (length(miss_hs) > 0) stop("HashSolo table missing columns: ", paste(miss_hs, collapse = ", "), " :: ", hs_path)

hs_tab <- hs %>%
  transmute(
    capture = as.character(HS___capture__),
    sample_short = as.character(HS_Sample_short),
    n_cells = as.numeric(n_cells),
    method = "HashSolo"
  ) %>%
  filter(!is.na(capture), capture != "", !is.na(sample_short), sample_short != "", is.finite(n_cells))

write_tsv(
  hs_tab %>% transmute(HS___capture__ = capture, HS_Sample_short = sample_short, n_cells),
  file.path(outdir, "hashsolo_counts_per_capture_sample.tsv")
)
msg("Wrote HashSolo table: ", file.path(outdir, "hashsolo_counts_per_capture_sample.tsv"))

# -------------------------
# 2) Build genetic demultiplexing table from *_clusters_with_patient.tsv
# -------------------------
gen_files <- list.files(gen_root, pattern="*_clusters_with_patient\\.tsv$", recursive=TRUE, full.names=TRUE)
if (length(gen_files) == 0) stop("No *_clusters_with_patient.tsv under: ", gen_root)

read_one_gen <- function(f) {
  cap <- basename(dirname(f))  # parent folder is capture (C1, C6, ...)
  df <- readr::read_tsv(f, show_col_types = FALSE, col_types = cols(.default = col_character()))
  req <- c("status", "patient_id_from_HTO")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Genetic file missing columns: ", paste(miss, collapse=", "), " :: ", f)

  df %>%
    filter(status == "singlet", !is.na(patient_id_from_HTO), patient_id_from_HTO != "") %>%
    count(gen_dem_capture = cap, gen_dem_Sample_short = patient_id_from_HTO, name = "n_cells")
}

gen_tab_wide <- purrr::map_dfr(sort(gen_files), read_one_gen) %>%
  group_by(gen_dem_capture, gen_dem_Sample_short) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  arrange(gen_dem_capture, gen_dem_Sample_short)

write_tsv(gen_tab_wide, file.path(outdir, "genetic_counts_per_capture_sample.tsv"))
msg("Wrote genetic table: ", file.path(outdir, "genetic_counts_per_capture_sample.tsv"))

gen_tab <- gen_tab_wide %>%
  transmute(
    capture = as.character(gen_dem_capture),
    sample_short = as.character(gen_dem_Sample_short),
    n_cells = as.numeric(n_cells),
    method = "Genetic"
  ) %>%
  filter(!is.na(capture), capture != "", !is.na(sample_short), sample_short != "", is.finite(n_cells))

# -------------------------
# 3) Combine + plot
# IMPORTANT FIX: do NOT 'complete()' across captures, otherwise every capture
# gets every sample (with zeros) and you see all samples in every facet.
# We only expand methods within each (capture, sample_short) that exists.
# -------------------------
comb_raw <- bind_rows(hs_tab, gen_tab) %>%
  group_by(capture, sample_short, method) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

# only keep samples that exist in that capture in at least one method
samples_by_cap <- comb_raw %>% distinct(capture, sample_short)

comb <- samples_by_cap %>%
  tidyr::expand_grid(method = c("HashSolo","Genetic")) %>%
  left_join(comb_raw, by = c("capture","sample_short","method")) %>%
  mutate(
    n_cells = dplyr::coalesce(n_cells, 0),
    method = factor(method, levels = c("HashSolo","Genetic"))
  )

write_tsv(comb, file.path(outdir, "combined_counts_long.tsv"))
msg("Wrote combined long table: ", file.path(outdir, "combined_counts_long.tsv"))

# Order samples within each capture by total (HashSolo+Genetic)
ord <- comb %>%
  group_by(capture, sample_short) %>%
  summarise(total = sum(n_cells), .groups = "drop") %>%
  arrange(capture, desc(total)) %>%
  mutate(sample_capture = paste(capture, sample_short, sep="__"))

levels_sc <- ord$sample_capture

plot_df <- comb %>%
  mutate(
    sample_capture = paste(capture, sample_short, sep="__"),
    sample_capture = factor(sample_capture, levels = levels_sc)
  )

p <- ggplot(plot_df, aes(x = sample_capture, y = n_cells, fill = method)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  facet_wrap(~ capture, scales = "free_x") +
  scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
  labs(
    title = "Cells per sample per capture: HashSolo vs Genetic demultiplexing",
    x = "Sample (short)",
    y = "Number of cells",
    fill = "Method"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

png_out <- file.path(outdir, "HS_vs_genetic_counts_per_capture.png")
pdf_out <- file.path(outdir, "HS_vs_genetic_counts_per_capture.pdf")
ggsave(png_out, p, width = 14, height = 6, dpi = 300)
ggsave(pdf_out, p, width = 14, height = 6)

msg("Wrote plot: ", png_out)
msg("Wrote plot: ", pdf_out)

totals <- comb %>%
  group_by(capture, method) %>%
  summarise(n_cells = sum(n_cells), .groups="drop") %>%
  arrange(capture, method)
msg("Totals per capture x method:")
print(totals, n = nrow(totals))
RS

log "DONE. Output folder: $OUTDIR"
log "Tables:"
log "  $OUTDIR/hashsolo_counts_per_capture_sample.tsv"
log "  $OUTDIR/genetic_counts_per_capture_sample.tsv"
log "  $OUTDIR/combined_counts_long.tsv"
log "Plots:"
log "  $OUTDIR/HS_vs_genetic_counts_per_capture.png"
log "  $OUTDIR/HS_vs_genetic_counts_per_capture.pdf"
