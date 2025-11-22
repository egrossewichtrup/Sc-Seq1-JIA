#!/usr/bin/env bash
#SBATCH --job-name=plot_souporcell_C1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --chdir=/hpc/dla_lti/yermanos_group/GRO11679
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=egrossewichtrup@umcutrecht.nl

set -euo pipefail

log(){ echo "[$(date +%F_%T)] $*"; }

########################
# USER SETTINGS       #
########################

# Capture to plot
CAPTURE="C1"

# Path to genetic demultiplexing run root
SP_ROOT="/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing/demux/genetic_demultiplexing_20251119_144453"

# Subfolder created by the annotate_souporcell job
ANNOT_SUBDIR="souporcell_with_patient_20251122_095429"

# R binary (same one as used in your other pipelines)
VENV_ROOT="/hpc/dla_lti/yermanos_group/software"
R_BIN="$VENV_ROOT/venvs4/r/bin/Rscript"

########################
# END USER SETTINGS    #
########################

log "Starting plotting job for capture ${CAPTURE}"

if [ ! -x "$R_BIN" ]; then
  log "FATAL: Rscript not found or not executable at $R_BIN"
  exit 1
fi

CAPTURE_DIR="${SP_ROOT}/${ANNOT_SUBDIR}/${CAPTURE}"
CLUSTERS_FILE="${CAPTURE_DIR}/${CAPTURE}_clusters_with_patient.tsv"

if [ ! -f "$CLUSTERS_FILE" ]; then
  log "FATAL: clusters_with_patient file not found:"
  log "       $CLUSTERS_FILE"
  exit 1
fi

log "Using clusters file: $CLUSTERS_FILE"

# TMPDIR for this job
PROJ_MAIN="/hpc/dla_lti/yermanos_group/GRO11679"
JOBTMP="${PROJ_MAIN}/tmp/plot_souporcell_${SLURM_JOB_ID:-noid}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$JOBTMP"
export TMPDIR="$JOBTMP"
export TMP="$JOBTMP"
export TEMP="$JOBTMP"
log "TMPDIR set to $JOBTMP"

# Write R script into project scripts dir
SCRIPT_DIR="${PROJ_MAIN}/scripts"
mkdir -p "$SCRIPT_DIR"
RSCRIPT_FILE="${SCRIPT_DIR}/plot_souporcell_barplot_${CAPTURE}.R"

log "Writing R script to $RSCRIPT_FILE"

cat > "$RSCRIPT_FILE" <<'RS'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# --- read env vars ----
capture      <- Sys.getenv("CAPTURE", "C1")
clusters_file <- Sys.getenv("CLUSTERS_FILE")
capture_dir   <- dirname(clusters_file)

message("[INFO] Capture:       ", capture)
message("[INFO] Clusters file: ", clusters_file)

if (clusters_file == "" || !file.exists(clusters_file)) {
  stop("clusters_with_patient.tsv not found at: ", clusters_file)
}

# --- load data ----
df <- read_tsv(clusters_file, show_col_types = FALSE)

needed <- c("barcode", "status", "assignment", "patient_id_from_HTO")
missing <- setdiff(needed, names(df))
if (length(missing) > 0) {
  stop("Missing columns in ", clusters_file, ": ",
       paste(missing, collapse = ", "))
}

# --- filter to singlets with mapped patient ----
df_filt <- df %>%
  filter(
    status == "singlet",
    !is.na(patient_id_from_HTO)
  )

if (nrow(df_filt) == 0) {
  stop("No singlet cells with patient_id_from_HTO found for capture ", capture)
}

message("[INFO] Cells used for plotting (singlets with patient ID): ",
        nrow(df_filt))

# --- summarise counts per cluster & patient ----
counts <- df_filt %>%
  count(assignment, patient_id_from_HTO, name = "n_cells") %>%
  group_by(assignment) %>%
  mutate(n_total_cluster = sum(n_cells)) %>%
  ungroup()

# order clusters numerically where possible
counts <- counts %>%
  mutate(
    assignment = factor(assignment,
                        levels = sort(unique(assignment)))
  )

# --- barplot: absolute number of cells ----
p_cells <- ggplot(counts,
                  aes(x = assignment,
                      y = n_cells,
                      fill = patient_id_from_HTO)) +
  geom_col() +
  labs(
    title = paste0("Souporcell genetic clusters for capture ", capture),
    x     = "Souporcell cluster (assignment)",
    y     = "Number of cells (singlets)",
    fill  = "Patient (HTO)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

out_png_cells <- file.path(
  capture_dir,
  paste0(capture, "_souporcell_barplot_cells_per_cluster_by_patient.png")
)

ggsave(out_png_cells, p_cells, width = 7, height = 4.5, dpi = 300)
message("[INFO] Saved absolute-count barplot to: ", out_png_cells)

# --- barplot: fraction per cluster ----
counts_frac <- counts %>%
  group_by(assignment) %>%
  mutate(frac = n_cells / n_total_cluster) %>%
  ungroup()

p_frac <- ggplot(counts_frac,
                 aes(x = assignment,
                     y = frac,
                     fill = patient_id_from_HTO)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    title = paste0("Souporcell cluster composition for capture ", capture),
    x     = "Souporcell cluster (assignment)",
    y     = "Fraction of cells per cluster",
    fill  = "Patient (HTO)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

out_png_frac <- file.path(
  capture_dir,
  paste0(capture, "_souporcell_barplot_fraction_per_cluster_by_patient.png")
)

ggsave(out_png_frac, p_frac, width = 7, height = 4.5, dpi = 300)
message("[INFO] Saved fraction barplot to: ", out_png_frac)

message("[DONE] Finished plotting for capture ", capture)
RS

chmod +x "$RSCRIPT_FILE"

export CAPTURE CLUSTERS_FILE

log "Running R script with:"
log "  CAPTURE       = $CAPTURE"
log "  CLUSTERS_FILE = $CLUSTERS_FILE"

"$R_BIN" "$RSCRIPT_FILE"

log "R script finished successfully"
log "Barplots are in: $CAPTURE_DIR"
log "Temporary directory was: $JOBTMP"
