#!/usr/bin/env bash
#SBATCH --job-name=DEMUX_PIPELINE
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=egrossewichtrup@umcutrecht.nl
#SBATCH --chdir=/hpc/dla_lti/yermanos_group/GRO11679

set -euo pipefail

# -------- logger early (needed before first log call) ----------
log(){ echo "[$(date +%F_%T)] $*"; }

#############################
# USER SETTINGS (edit here) #
#############################
# Project root for THIS run (all outputs land here)
PROJ="/hpc/dla_lti/yermanos_group/GRO11679/demultiplexing"
mkdir -p "$PROJ"
cd "$PROJ"

# Existing Miniforge/Python/R envs live here (do NOT put in PROJ unless you copied them)
VENV_ROOT="/hpc/dla_lti/yermanos_group/software"

# Cell Ranger multi ROOT containing your dated folder with mixed names & spaces
CR_ROOT="/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/cellranger_multi_2025-10-30"

# Interpreters (from VENV_ROOT, not PROJ)
MF3="$VENV_ROOT/venvs4/miniforge3"
PYENV="$VENV_ROOT/venvs4/py"
R_BIN="$VENV_ROOT/venvs4/r/bin/Rscript"

# HashSolo defaults
HASH_SOLO_PRIORS="0.01,0.80,0.19"
SINGLET_CONFIDENCE_MIN="0.90"

#########################################
# CAPTURE CONFIG (edit values below)    #
# Columns (tab-separated):
# cap   runfolder                                   htos                 map_short                                  map_full
# - cap: your canonical ID used in labels/plots
# - runfolder: exact per-sample-outs folder name ("multi_*") as produced by Cell Ranger
# - htos: comma list of HTO feature names as they appear in the Antibody Capture matrix
# - map_short: semicolon list of HTO=ShortSample
# - map_full:  semicolon list of HTO=FullSample
#########################################

# 1) Write config to a temporary location under $PROJ (RUNROOT not known yet)
TS="$(date +%Y%m%d_%H%M%S)"
TMP_CFG_DIR="$PROJ/tmp_cfg_${TS}"
mkdir -p "$TMP_CFG_DIR"
TMP_CAPTURE_CFG_TSV="$TMP_CFG_DIR/captures.tsv"

cat > "$TMP_CAPTURE_CFG_TSV" <<'CFG'
C1	multi_C1_ScSeq1-Capture1-HC001-HC002-10-07-2025	HTO1,HTO2	HTO1=HC001;HTO2=HC002	HTO1=HC001_Sc-Seq1_PBMC_Bcells_10.07.2025;HTO2=HC002_Sc-Seq1_PBMC_Bcells_10.07.2025
C6	multi_C6_ScSeq1-Capture6-P003-HC003-21-07-2025	HTO1,HTO2	HTO1=P003;HTO2=HC003	HTO1=P003_Sc-Seq1_PBMC_Bcells_21.07.2025;HTO2=HC003_Sc-Seq1_PBMC_Bcells_21.07.2025
C8	multi_C8_ScSeq1-Capture8-P005-P006-24-07-2025	HTO3,HTO4	HTO3=P005;HTO4=P006	HTO3=P005_Sc-Seq1_PBMC_Bcells_24.07.2025;HTO4=P006_Sc-Seq1_PBMC_Bcells_24.07.2025
C9	multi_C9_ScSeq1-Capture9-P002blood-P007-P008-30-07-2025	HTO1,HTO2,HTO3	HTO1=P002;HTO2=P007;HTO3=P008	HTO1=P002_Sc-Seq1_PPBMC_Bcells_30.07.2025;HTO2=P007_Sc-Seq1_PBMC_Bcells_30.07.2025;HTO3=P008_Sc-Seq1_PBMC_Bcells_30.07.2025
C10	multi_C10_ScSeq1-Capture10-P009-P010-31-07-2025	HTO3,HTO4	HTO3=P009;HTO4=P010	HTO3=P009_Sc-Seq1_PBMC_Bcells_31.07.2025;HTO4=P010_Sc-Seq1_PBMC_Bcells_31.07.2025
CFG
log "Wrote temporary capture config: $TMP_CAPTURE_CFG_TSV"

# 2) Derive capture IDs and tag from the temporary config
CAPTURE_IDS="$(awk -F'\t' 'NF>=2{print $1}' "$TMP_CAPTURE_CFG_TSV" | paste -sd' ' -)"
TAG_JOINED="$(echo "$CAPTURE_IDS" | sed 's/ /_/g')"

########################
# Derived paths
########################
RUNROOT="$PROJ/demux/${TS}_demux_${TAG_JOINED}"
LOGS="$RUNROOT/logs"
SCRIPTS="$RUNROOT/scripts"
RSCRIPTS="$SCRIPTS/r"
PYSCRIPTS="$SCRIPTS/py"
SBATCHS="$SCRIPTS/sbatch"
OUT_TOP="$RUNROOT/out"
CFGDIR="$RUNROOT/config"
mkdir -p "$RUNROOT" "$LOGS" "$RSCRIPTS" "$PYSCRIPTS" "$SBATCHS" "$OUT_TOP" "$CFGDIR"

# 3) Move the config into the final run folder and re-point the variable
CAPTURE_CFG_TSV="$CFGDIR/captures.tsv"
mv "$TMP_CAPTURE_CFG_TSV" "$CAPTURE_CFG_TSV"
log "Capture config now at: $CAPTURE_CFG_TSV"

# Seurat / HashSolo / QC output roots
SEUROOT="$PROJ/seurat/${TS}_seurat_merge_${TAG_JOINED}"
HASHROOT="$PROJ/hashsolo/${TS}_hashsolo_${TAG_JOINED}"
QCCOUNTS="$PROJ/qc/step4_1_counts_${TS}"
FILTERROOT="$PROJ/qc/step4_2_filter_${TS}"
SPLITROOT="$PROJ/seurat/${TS}_step4_3_split"
VDJROOT="$PROJ/vdj/${TS}_VDJ_PER_SAMPLE"

log "DEMUX driver starting. Run root: $RUNROOT"

##################################
# Preflight
##################################
# Rscript finder (use venv, then PATH)
if [ -x "$R_BIN" ]; then :; elif command -v Rscript >/dev/null 2>&1; then R_BIN="$(command -v Rscript)"; else
  log "ERROR: No Rscript found. Provide at $R_BIN or load it via modules."; exit 2; fi

# Check Miniforge & Python env
if [ ! -f "$MF3/etc/profile.d/conda.sh" ]; then log "ERROR: Miniforge conda.sh not found at $MF3/etc/profile.d/conda.sh"; exit 2; fi
if [ ! -d "$PYENV" ]; then log "ERROR: Python env not found at $PYENV"; exit 2; fi

##################################
# Build H5/VDJ indices (resolves paths even with spaces)
##################################
H5_INDEX_TSV="$CFGDIR/h5_index.tsv"
VDJ_INDEX_TSV="$CFGDIR/vdj_index.tsv"
echo -e "cap\trunfolder\th5_path" > "$H5_INDEX_TSV"
echo -e "cap\trunfolder\tvdj_dir\tcsv\tfasta" > "$VDJ_INDEX_TSV"

while IFS=$'\t' read -r CAP RUNFOLDER HTOS MAPS MAPF; do
  [ -z "${CAP:-}" ] && continue
  FOUND="$(find "$CR_ROOT" -type d -name "$RUNFOLDER" -print -quit || true)"
  if [ -z "$FOUND" ]; then
    log "MISSING runfolder: $RUNFOLDER (cap $CAP)"; continue
  fi
  OUTS="$FOUND/outs/per_sample_outs/$RUNFOLDER"
  H5="$OUTS/count/sample_filtered_feature_bc_matrix.h5"
  if [ -s "$H5" ]; then
    echo -e "${CAP}\t${RUNFOLDER}\t${H5}" >> "$H5_INDEX_TSV"
    log "OK H5 ($CAP): $H5"
  else
    log "MISSING H5 ($CAP): $H5"
  fi
  VDJ_DIR="$OUTS/vdj_b"
  CSV="$VDJ_DIR/filtered_contig_annotations.csv"
  FAS="$VDJ_DIR/filtered_contig.fasta"
  if [ -d "$VDJ_DIR" ] && [ -s "$CSV" ] && [ -s "$FAS" ]; then
    echo -e "${CAP}\t${RUNFOLDER}\t${VDJ_DIR}\t${CSV}\t${FAS}" >> "$VDJ_INDEX_TSV"
  else
    log "[WARN] VDJ missing or incomplete for $CAP at $VDJ_DIR"
  fi
done < "$CAPTURE_CFG_TSV"

# Abort early if any H5 missing
if [ "$(awk 'END{print NR}' "$H5_INDEX_TSV")" -le 1 ]; then
  log "ERROR: No valid H5s discovered. Check CR_ROOT and runfolder names."; exit 3
fi

##########################################
# STEP 1 — Write + run Seurat merge (R)  #
##########################################
LOAD_R="$RSCRIPTS/load_merge_seurat.R"
cat > "$LOAD_R" <<'RS'
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  libs <- c("Seurat","Matrix","hdf5r","readr","dplyr")
  missing <- setdiff(libs, rownames(installed.packages()))
  if (length(missing) > 0) stop("Missing R packages: ", paste(missing, collapse=", "))
  library(Seurat); library(Matrix); library(hdf5r); library(readr); library(dplyr)
})
CR_ROOT <- Sys.getenv("CR_ROOT")
SEUROOT <- Sys.getenv("SEUROOT")
H5INDEX <- Sys.getenv("H5INDEX")
SEU_TAG <- Sys.getenv("SEU_TAG")

if (any(nchar(c(CR_ROOT,SEUROOT,H5INDEX))==0)) stop("Missing env: CR_ROOT/SEUROOT/H5INDEX")
dir.create(SEUROOT, recursive=TRUE, showWarnings=FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), paste0(...)))

h5tab <- readr::read_tsv(H5INDEX, show_col_types=FALSE, col_types=cols(.default=col_character()))
stopifnot(all(c("cap","runfolder","h5_path") %in% names(h5tab)))
if (nrow(h5tab)==0) stop("No H5 entries present.")

seurat_list <- list(); qc_lines <- c()
for (i in seq_len(nrow(h5tab))) {
  cap <- h5tab$cap[i]; h5 <- h5tab$h5_path[i]
  if (!file.exists(h5)) { warning("H5 not found: ", h5); next }
  msg("Loading ", h5)
  data <- Read10X_h5(h5)
  stopifnot("Gene Expression" %in% names(data))
  seu <- CreateSeuratObject(counts=data$`Gene Expression`, project=paste0("multi_",cap))
  if ("Antibody Capture" %in% names(data)) {
    ab <- data$`Antibody Capture`
    seu[["HTO"]] <- CreateAssayObject(counts = ab)
    seu <- NormalizeData(seu, assay="HTO", normalization.method="CLR", margin=2, verbose=FALSE)
    qc_lines <- c(qc_lines, sprintf("%s: cells=%d HTO_features=%d", cap, ncol(seu), nrow(ab)))
  } else qc_lines <- c(qc_lines, sprintf("%s: WARNING no Antibody Capture", cap))
  seurat_list[[paste0("multi_",cap)]] <- seu
}
if (length(seurat_list)==0) stop("No Seurat objects created.")
combined <- if (length(seurat_list)==1) seurat_list[[1]] else Reduce(function(x,y) merge(x,y, add.cell.ids=c(x@project.name, y@project.name)), seurat_list)
tag <- if (nchar(SEU_TAG)>0) SEU_TAG else "MERGED"
rds_path <- file.path(SEUROOT, paste0("combined_multi_seurat_clean_", tag, ".rds"))
saveRDS(combined, rds_path)
msg("Saved combined Seurat: ", rds_path)

if ("HTO" %in% Assays(combined)) {
  hto_raw <- GetAssayData(combined, assay="HTO", slot="counts")
  csv_path <- file.path(SEUROOT, paste0("hto_counts_raw_", tag, ".csv"))
  write.csv(as.matrix(hto_raw), csv_path, quote=FALSE)
  msg("Exported HTO raw CSV: ", csv_path)
} else msg("WARNING: No HTO assay present; no CSV exported.")
dir.create(file.path(SEUROOT,"logs"), showWarnings=FALSE, recursive=TRUE)
writeLines(qc_lines, con=file.path(SEUROOT,"logs",paste0("qc_summary_",tag,".txt")))
RS
chmod +x "$LOAD_R"

export CR_ROOT SEUROOT
export H5INDEX="$H5_INDEX_TSV"
export SEU_TAG="${TS}_${TAG_JOINED}"
log "STEP 1: Seurat merge"
"$R_BIN" "$LOAD_R" 2>&1 | tee "$LOGS/step1_seurat_merge.log"

# Discover produced files
MERGED_RDS="$(ls -1 "$SEUROOT"/combined_multi_seurat_clean_*_${TAG_JOINED}.rds 2>/dev/null | tail -n1 || true)"
[ -z "$MERGED_RDS" ] && MERGED_RDS="$(ls -1 "$SEUROOT"/combined_multi_seurat_clean_*_MERGED.rds 2>/dev/null | tail -n1 || true)"
HTO_CSV="$(ls -1 "$SEUROOT"/hto_counts_raw_*_${TAG_JOINED}.csv 2>/dev/null | tail -n1 || true)"
[ -z "$HTO_CSV" ] && HTO_CSV="$(ls -1 "$SEUROOT"/hto_counts_raw_*_MERGED.csv 2>/dev/null | tail -n1 || true)"
[ -s "$MERGED_RDS" ] || { log "FATAL: merged RDS not found"; exit 10; }
[ -s "$HTO_CSV" ]    || { log "FATAL: HTO CSV not found"; exit 11; }
log "Found RDS: $MERGED_RDS"
log "Found CSV: $HTO_CSV"

########################################################
# STEP 2 — HashSolo: CSV→H5AD + per-capture demux (Py) #
########################################################
HS_PY="$PYSCRIPTS/hashsolo_run.py"
cat > "$HS_PY" <<'PY'
import os, re, itertools, numpy as np, pandas as pd
import anndata as ad, scanpy as sc
import scanpy.external as sce
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

csv = os.environ["HTO_CSV"]
outd = os.environ["HASHROOT"]
cfg_tsv = os.environ["CAPTURE_CFG_TSV"]
priors = tuple(float(x) for x in os.environ.get("HASH_SOLO_PRIORS","0.01,0.80,0.19").split(","))
singlet_min = float(os.environ.get("SINGLET_CONFIDENCE_MIN","0.90"))

os.makedirs(outd, exist_ok=True)
print(f"[INFO] Reading HTO CSV: {csv}")
df = pd.read_csv(csv, index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
ad_all = ad.AnnData(df.T.astype(np.float64))
ad_all.obs_names_make_unique(); ad_all.var_names_make_unique()
ad_all.write(os.path.join(outd, "hto_counts_raw.h5ad"))
print(f"[INFO] Wrote {os.path.join(outd,'hto_counts_raw.h5ad')} shape={ad_all.shape}")

cfg = pd.read_csv(cfg_tsv, sep="\t", header=None, names=["cap","runfolder","htos","map_short","map_full"])
cfg["cap"] = cfg["cap"].astype(str)

def cells_of_cap(names, cap):
    pat = re.compile(rf"(?:^|_)multi_{re.escape(cap)}_")
    return [c for c in names if pat.search(c)]

obs_parts = []
for _, row in cfg.iterrows():
    cap = row["cap"]
    keep = cells_of_cap(ad_all.obs_names, cap)
    print(f"[INFO] {cap}: cells={len(keep)}")
    if not keep:
        print(f"[WARN] {cap}: no cells; skip"); continue
    sub = ad_all[keep,:].copy()

    want = [h.strip() for h in str(row["htos"]).split(",") if h.strip() != ""]
    present = [h for h in want if h in sub.var_names]
    missing = [h for h in want if h not in present]
    if missing:
        print(f"[WARN] {cap}: missing HTO features not in matrix: {missing}")
    if len(present) < 2:
        print(f"[WARN] {cap}: <2 HTOs present ({present}); skip"); continue

    for h in present:
        j = list(sub.var_names).index(h)
        col = sub.X[:, j]
        try: v = col.toarray().ravel().astype(float)
        except AttributeError: v = np.asarray(col).ravel().astype(float)
        sub.obs[h] = v

    kwargs = {"number_of_noise_barcodes": 1} if len(present)==2 else {}
    print(f"[INFO] {cap}: running HashSolo HTOs={present} priors={priors} extra={kwargs}")
    sce.pp.hashsolo(sub, present, priors=priors, **kwargs)

    cap_dir = os.path.join(outd, cap); os.makedirs(cap_dir, exist_ok=True)
    sub.write(os.path.join(cap_dir, f"hashsoloed_{cap}.h5ad"))

    pdf = os.path.join(cap_dir, f"hashsolo_qc_{cap}.pdf")
    with PdfPages(pdf) as pp:
        vc = sub.obs["Classification"].value_counts().sort_index()
        plt.figure(figsize=(6,4)); vc.plot(kind="bar"); plt.title(f"{cap} — HashSolo classes")
        plt.ylabel("Cells"); plt.xlabel("Class"); plt.tight_layout(); pp.savefig(); plt.close()

        if "singlet_hypothesis_probability" in sub.obs.columns:
            plt.figure(figsize=(6,4))
            sub.obs["singlet_hypothesis_probability"].plot(kind="hist", bins=50)
            plt.title(f"{cap} — Singlet prob"); plt.xlabel("Probability"); plt.ylabel("Cells")
            plt.tight_layout(); pp.savefig(); plt.close()

        htos = [h for h in present if h in sub.obs.columns]
        for x,y in itertools.combinations(htos, 2):
            plt.figure(figsize=(5,5))
            X = np.log10(sub.obs[x].values + 1); Y = np.log10(sub.obs[y].values + 1)
            plt.scatter(X,Y, s=2, alpha=0.4)
            plt.xlabel(f"log10({x}+1)"); plt.ylabel(f"log10({y}+1)"); plt.title(f"{cap} — {x} vs {y}")
            plt.tight_layout(); pp.savefig(); plt.close()

    obs = sub.obs.copy(); obs["__capture__"] = cap; obs_parts.append(obs)

if not obs_parts:
    raise SystemExit("No HashSolo obs collected.")
full = pd.concat(obs_parts, axis=0)
full.to_csv(os.path.join(outd, "hashsolo_cells_full.csv"), index=True, index_label="cell_id")
print(f"[OK] hashsolo_cells_full.csv rows={full.shape[0]}")
PY
chmod +x "$HS_PY"

log "STEP 2: HashSolo (CSV→H5AD + per-capture run)"
source "$MF3/etc/profile.d/conda.sh"
conda activate "$PYENV"
export HTO_CSV HASHROOT
export CAPTURE_CFG_TSV="$CAPTURE_CFG_TSV"
export HASH_SOLO_PRIORS SINGLET_CONFIDENCE_MIN
python "$HS_PY" 2>&1 | tee "$LOGS/step2_hashsolo.log"

HS_CSV="$HASHROOT/hashsolo_cells_full.csv"
[ -s "$HS_CSV" ] || { log "FATAL: HashSolo CSV missing"; exit 20; }

###########################################################
# STEP 3 — Attach HS_* back to Seurat & map sample names  #
###########################################################
ATTACH_R="$RSCRIPTS/attach_hashsolo.R"
cat > "$ATTACH_R" <<'RS'
suppressPackageStartupMessages({ library(Seurat); library(readr); library(dplyr); library(tibble); library(stringr) })
seurat_rds <- Sys.getenv("MERGED_RDS")
hs_csv     <- Sys.getenv("HS_CSV")
stopifnot(file.exists(seurat_rds), file.exists(hs_csv))
so <- readRDS(seurat_rds)
hs <- read_csv(hs_csv, show_col_types=FALSE)
if ("...1" %in% names(hs) && !"cell_id" %in% names(hs)) hs <- rename(hs, cell_id = `...1`)
stopifnot("cell_id" %in% names(hs))
keep_cols <- setdiff(names(hs), "cell_id")
hs <- hs %>% rename_with(~ paste0("HS_", .), all_of(keep_cols)) %>%
  mutate(HS_is_singlet = HS_Classification %in% paste0("HTO", 1:128))
hs_df <- hs %>% column_to_rownames("cell_id") %>% as.data.frame()
hs_df <- hs_df[colnames(so), , drop=FALSE]
so <- AddMetaData(so, metadata=hs_df)
out_rds <- sub("\\.rds$", "_with_hashsolo.rds", seurat_rds)
saveRDS(so, out_rds); cat("[OK] Saved: ", out_rds, "\n", sep="")
RS
chmod +x "$ATTACH_R"

export MERGED_RDS HS_CSV
log "STEP 3a: Attach HS_* to Seurat"
"$R_BIN" "$ATTACH_R" 2>&1 | tee "$LOGS/step3a_attach.log"
WITH_HS_RDS="${MERGED_RDS%.rds}_with_hashsolo.rds"
[ -s "$WITH_HS_RDS" ] || { log "FATAL: _with_hashsolo.rds not found"; exit 30; }

# Map names from config TSV
MAPNAMES_R="$RSCRIPTS/map_names.R"
cat > "$MAPNAMES_R" <<'RS'
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(readr); library(stringr) })
in_rds <- Sys.getenv("IN_RDS")
cfg_tsv <- Sys.getenv("CAPTURE_CFG_TSV")
so <- readRDS(in_rds); stopifnot(inherits(so,"Seurat"))
req <- c("HS_Classification","HS___capture__"); stopifnot(all(req %in% colnames(so@meta.data)))

cfg <- readr::read_tsv(cfg_tsv, col_names=c("cap","runfolder","htos","map_short","map_full"), show_col_types=FALSE)

parse_map <- function(s){
  if(is.na(s) || s=="") return(setNames(character(0), character(0)))
  parts <- strsplit(s, ";")[[1]]
  kv <- strsplit(parts, "=")
  keys <- vapply(kv, function(x) trimws(x[1]), character(1))
  vals <- vapply(kv, function(x) trimws(x[2]), character(1))
  setNames(vals, keys)
}

map_short <- setNames(vector("list", nrow(cfg)), cfg$cap)
map_full  <- setNames(vector("list", nrow(cfg)), cfg$cap)
for (i in seq_len(nrow(cfg))) {
  map_short[[ cfg$cap[i] ]] <- parse_map(cfg$map_short[i])
  map_full [[ cfg$cap[i] ]] <- parse_map(cfg$map_full[i])
}

is_singlet <- so$HS_Classification %in% paste0("HTO", 1:128)
lookup <- function(cap, call, MAP){
  if(is.na(cap) || is.na(call)) return(NA_character_)
  if(!cap %in% names(MAP)) return(NA_character_)
  m <- MAP[[cap]]
  if(call %in% names(m)) m[[call]] else NA_character_
}

so$HS_Sample_full  <- ifelse(is_singlet, mapply(lookup, so$HS___capture__, so$HS_Classification, MoreArgs=list(MAP=map_full)),  NA_character_)
so$HS_Sample_short <- ifelse(is_singlet, mapply(lookup, so$HS___capture__, so$HS_Classification, MoreArgs=list(MAP=map_short)), NA_character_)

if ("HS_singlet_hypothesis_probability" %in% colnames(so@meta.data)) {
  so$HS_is_confident <- is_singlet & (so$HS_singlet_hypothesis_probability >= as.numeric(Sys.getenv("SINGLET_CONFIDENCE_MIN","0.90")))
} else {
  so$HS_is_confident <- is_singlet
}

out_rds <- sub("\\.rds$", "_named.rds", in_rds)
saveRDS(so, out_rds)
cat("[OK] Saved: ", out_rds, "\n", sep="")
RS
chmod +x "$MAPNAMES_R"

export IN_RDS="$WITH_HS_RDS" CAPTURE_CFG_TSV SINGLET_CONFIDENCE_MIN
log "STEP 3b: Map sample names"
"$R_BIN" "$MAPNAMES_R" 2>&1 | tee "$LOGS/step3b_mapnames.log"

NAMED_RDS="${WITH_HS_RDS%.rds}_named.rds"
[ -s "$NAMED_RDS" ] || { log "FATAL: named RDS not found"; exit 31; }

##########################################################
# STEP 4 — QC stacked counts (png) + STEP 5 filter 4.2  #
##########################################################
QC_R="$RSCRIPTS/qc_counts_png.R"
cat > "$QC_R" <<'RS'
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(readr); library(ggplot2); library(tidyr) })
in_rds <- Sys.getenv("IN_RDS"); outdir <- Sys.getenv("OUTDIR"); CONF_T <- as.numeric(Sys.getenv("CONF_T","0.90"))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
so <- readRDS(in_rds); md <- so@meta.data
need <- c("HS___capture__","HS_Classification","HS_Sample_short","HS_singlet_hypothesis_probability")
miss <- setdiff(need, colnames(md)); if (length(miss)>0) stop("Missing cols: ", paste(miss, collapse=", "))
md$HS_is_confident_effective <- (md$HS_Classification %in% c("Singlet","SINGLET","singlet", paste0("HTO",1:128))) &
                                is.finite(md$HS_singlet_hypothesis_probability) &
                                md$HS_singlet_hypothesis_probability >= CONF_T
if (sum(md$HS_is_confident_effective,na.rm=TRUE)==0) {
  md$HS_is_confident_effective <- !(md$HS_Classification %in% c("Doublet","doublet","Negative","negative"))
}
md_sel <- md %>% filter(HS_is_confident_effective %in% TRUE)
per_cap_sample <- md_sel %>%
  count(HS___capture__, HS_Sample_short, name="n_cells") %>%
  group_by(HS___capture__) %>%
  mutate(n_total_capture=sum(n_cells), frac_capture=n_cells/n_total_capture) %>% ungroup()
readr::write_csv(per_cap_sample, file.path(outdir, "counts_per_capture_sample.csv"))
readr::write_csv(per_cap_sample %>% count(HS___capture__, wt=n_cells, name="n_cells"), file.path(outdir, "counts_per_capture.csv"))
sample_order <- per_cap_sample %>% group_by(HS_Sample_short) %>% summarise(n=sum(n_cells), .groups="drop") %>% arrange(desc(n)) %>% pull(HS_Sample_short)
plot_df <- per_cap_sample %>% mutate(HS_Sample_short=factor(HS_Sample_short, levels=sample_order), HS___capture__=factor(HS___capture__))
p_counts <- ggplot(plot_df, aes(x=HS___capture__, y=n_cells, fill=HS_Sample_short)) + geom_col() +
  labs(x="Capture", y="Cells (confident)", fill="Sample", title=paste0("Cells per capture (prob ≥ ", CONF_T, ")")) + theme_bw() + theme(plot.title=element_text(hjust=.5))
p_frac <- ggplot(plot_df, aes(x=HS___capture__, y=frac_capture, fill=HS_Sample_short)) + geom_col(position="fill") +
  scale_y_continuous(labels=function(x) paste0(round(100*x),"%")) +
  labs(x="Capture", y="Fraction", fill="Sample", title=paste0("Fractions per capture (prob ≥ ", CONF_T, ")")) + theme_bw() + theme(plot.title=element_text(hjust=.5))
ggsave(file.path(outdir,"stacked_counts_per_capture.png"), p_counts, width=8, height=5, dpi=220)
ggsave(file.path(outdir,"stacked_fractions_per_capture.png"), p_frac, width=8, height=5, dpi=220)
RS
chmod +x "$QC_R"

export IN_RDS="$NAMED_RDS" OUTDIR="$QCCOUNTS/out" CONF_T="$SINGLET_CONFIDENCE_MIN"
log "STEP 4: QC stacked counts (PNGs)"
"$R_BIN" "$QC_R" 2>&1 | tee "$LOGS/step4_qc_counts.log"

FILTER_R="$RSCRIPTS/filter_confident_4_2.R"
cat > "$FILTER_R" <<'RS'
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(readr) })
in_rds <- Sys.getenv("IN_RDS"); outroot <- Sys.getenv("OUTROOT")
dir.create(outroot, recursive=TRUE, showWarnings=FALSE)
so <- readRDS(in_rds); md <- so[[]]
need <- c("HS_Classification","HS_is_confident","HS___capture__","HS_Sample_short")
miss <- setdiff(need, colnames(md)); if (length(miss)>0) stop("Missing: ", paste(miss,collapse=","))
is_doublet  <- md$HS_Classification %in% c("Doublet","doublet")
is_negative <- md$HS_Classification %in% c("Negative","negative")
is_conf     <- md$HS_is_confident %in% c(TRUE,"TRUE","True")
keep <- is_conf & !is_doublet & !is_negative
so_f <- subset(so, cells=colnames(so)[keep])
out_rds <- file.path(outroot, "seurat_confident_singlets.rds"); saveRDS(so_f, out_rds)
readr::write_csv(md |> mutate(keep=keep) |> count(HS___capture__, keep, name="n"), file.path(outroot, "counts_by_capture_flag.csv"))
cat("[OK] Saved: ", out_rds, "\n", sep="")
RS
chmod +x "$FILTER_R"

export IN_RDS="$NAMED_RDS" OUTROOT="$FILTERROOT"
log "STEP 5: Filter confident singlets (4.2)"
"$R_BIN" "$FILTER_R" 2>&1 | tee "$LOGS/step5_filter_4_2.log"
CONF_RDS="$FILTERROOT/seurat_confident_singlets.rds"
[ -s "$CONF_RDS" ] || { log "FATAL: confident singlets RDS missing"; exit 40; }

#####################################################
# STEP 6 — Split per capture + per sample (4.3)    #
#####################################################
SPLIT_R="$RSCRIPTS/split_4_3.R"
cat > "$SPLIT_R" <<'RS'
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(readr) })
input_rds <- Sys.getenv("INPUT_RDS"); outroot <- Sys.getenv("OUTROOT")
dir.create(outroot, recursive=TRUE, showWarnings=FALSE)
so <- readRDS(input_rds)
need <- c("HS___capture__","HS_Sample_short","HS_is_confident"); miss <- setdiff(need, colnames(so@meta.data))
if (length(miss)>0) stop("Missing: ", paste(miss,collapse=","))
if (!all(so@meta.data$HS_is_confident)) stop("Object must contain confident singlets only.")
out_seurat <- file.path(outroot, "seurat"); out_barcodes <- file.path(outroot, "barcodes"); out_idx <- file.path(outroot, "indices")
dir.create(out_seurat, TRUE, TRUE); dir.create(out_barcodes, TRUE, TRUE); dir.create(out_idx, TRUE, TRUE)
by_cap <- SplitObject(so, split.by="HS___capture__")
idx_cap <- list(); idx_samp <- list()
for (cap in sort(names(by_cap))) {
  so_cap <- by_cap[[cap]]; cap_dir <- file.path(out_seurat, cap); dir.create(cap_dir, TRUE, TRUE)
  cap_rds <- file.path(cap_dir, sprintf("seurat_capture_%s.rds", cap)); saveRDS(so_cap, cap_rds, compress="xz")
  bc_dir <- file.path(out_barcodes, "by_capture"); dir.create(bc_dir, TRUE, TRUE)
  cap_bc <- file.path(bc_dir, sprintf("%s_barcodes.csv", cap)); writeLines(colnames(so_cap), cap_bc)
  idx_cap[[cap]] <- data.frame(capture=cap, n_cells=ncol(so_cap), rds_path=cap_rds, barcodes_csv=cap_bc)
  by_s <- SplitObject(so_cap, split.by="HS_Sample_short")
  s_dir <- file.path(out_seurat, "by_sample", cap); s_bc_dir <- file.path(out_barcodes, "by_sample", cap)
  dir.create(s_dir, TRUE, TRUE); dir.create(s_bc_dir, TRUE, TRUE)
  for (s in sort(names(by_s))) {
    so_s <- by_s[[s]]; s_sanit <- gsub("[^A-Za-z0-9._-]+","_", s)
    s_rds <- file.path(s_dir, sprintf("seurat_%s_%s.rds", cap, s_sanit)); saveRDS(so_s, s_rds, compress="xz")
    s_bc  <- file.path(s_bc_dir, sprintf("%s_%s_barcodes.csv", cap, s_sanit)); writeLines(colnames(so_s), s_bc)
    idx_samp[[paste(cap,s)]] <- data.frame(capture=cap, sample=s, n_cells=ncol(so_s), rds_path=s_rds, barcodes_csv=s_bc)
  }
}
readr::write_csv(dplyr::bind_rows(idx_cap), file.path(out_idx, "index_per_capture.csv"))
readr::write_csv(dplyr::bind_rows(idx_samp), file.path(out_idx, "index_per_sample.csv"))
readr::write_csv(so@meta.data |> dplyr::count(HS___capture__, name="n_cells") |> arrange(HS___capture__), file.path(out_idx, "counts_per_capture.csv"))
readr::write_csv(so@meta.data |> dplyr::count(HS___capture__, HS_Sample_short, name="n_cells") |> arrange(HS___capture__, HS_Sample_short), file.path(out_idx, "counts_per_capture_sample.csv"))
cat("[OK] Split outputs at: ", outroot, "\n", sep="")
RS
chmod +x "$SPLIT_R"

export INPUT_RDS="$CONF_RDS" OUTROOT="$SPLITROOT"
log "STEP 6: Split per capture + sample"
"$R_BIN" "$SPLIT_R" 2>&1 | tee "$LOGS/step6_split_4_3.log"

###########################################################
# STEP 7 — VDJ subset per sample (sequential, lightweight)
###########################################################
MANIFEST="$VDJROOT/indices/step7_manifest.csv"
mkdir -p "$VDJROOT"/{indices,per_sample,logs}
echo "capture,sample,barcodes_csv,vdj_csv,vdj_fasta,out_dir" > "$MANIFEST"

while IFS=$'\t' read -r CAP RUNFOLDER VDJ_DIR CSV FAS; do
  [ "$CAP" = "cap" ] && continue
  BSDIR="$SPLITROOT/barcodes/by_sample/${CAP}"
  if [[ ! -s "$CSV" || ! -s "$FAS" || ! -d "$BSDIR" ]]; then
    log "[WARN] Skipping ${CAP} (missing VDJ or barcodes dir)"; continue
  fi
  find "$BSDIR" -maxdepth 1 -type f -name "${CAP}_*_barcodes.csv" | sort | while read -r BLIST; do
    BFILE="$(basename "$BLIST")"
    SAMPLE="${BFILE#${CAP}_}"; SAMPLE="${SAMPLE%_barcodes.csv}"
    ODIR="$VDJROOT/per_sample/${CAP}_${SAMPLE}"; mkdir -p "$ODIR"
    printf "%s,%s,%s,%s,%s,%s\n" "$CAP" "$SAMPLE" "$BLIST" "$CSV" "$FAS" "$ODIR" >> "$MANIFEST"
  done
done < <(tail -n +2 "$VDJ_INDEX_TSV")

log "VDJ manifest: $MANIFEST"
NTASKS=$(( $(wc -l < "$MANIFEST") - 1 ))
log "VDJ samples to subset: $NTASKS"

if [ "$NTASKS" -gt 0 ]; then
  tail -n +2 "$MANIFEST" | nl -w1 -s',' | while IFS=',' read -r IDX CAP SAMPLE BLIST CSV FAS ODIR; do
    log "[VDJ] [$IDX/$NTASKS] $CAP / $SAMPLE"
    HDR_LINE="$(head -n1 "$CSV")"
    BC_IDX=$(awk -F, -v H="$HDR_LINE" 'BEGIN{n=split(H,a,","); for(i=1;i<=n;i++){gsub(/^[[:space:]]+|[[:space:]]+$/,"",a[i]); if(a[i]=="barcode"){print i;exit}} }')
    CID_IDX=$(awk -F, -v H="$HDR_LINE" 'BEGIN{n=split(H,a,","); for(i=1;i<=n;i++){gsub(/^[[:space:]]+|[[:space:]]+$/,"",a[i]); if(a[i]=="contig_id"){print i;exit}} }')
    [ -n "${BC_IDX:-}" ] && [ -n "${CID_IDX:-}" ] || { log "[VDJ][ERROR] header columns missing"; exit 61; }

    BL_CORE="$ODIR/._bl_core.tmp"; CSV_CORE="$ODIR/._csv_core.tmp"; CID_TMP="$ODIR/._contigs.tmp"
    awk -F'_' 'NF{t=$NF; gsub(/\r$/,"",t); sub(/-1$/,"",t); if(t!="") print t}' "$BLIST" | sort -u > "$BL_CORE"
    awk -F, -v b="$BC_IDX" 'NR>1 && $b!=""{t=$b; sub(/-1$/,"",t); print t}' "$CSV" | sort -u > "$CSV_CORE"
    CSV_OUT="$ODIR/filtered_contig_annotations.${CAP}_${SAMPLE}.csv"
    FAS_OUT="$ODIR/filtered_contig.${CAP}_${SAMPLE}.fasta"
    awk -F, -v OFS=',' -v b="$BC_IDX" 'NR==FNR{keep[$1]=1;next} FNR==1{print;next}{t=$b; sub(/-1$/,"",t); if(t in keep) print}' "$BL_CORE" "$CSV" > "$CSV_OUT"
    awk -F, -v c="$CID_IDX" 'NR>1&&$c!=""{print $c}' "$CSV_OUT" | sort -u > "$CID_TMP"
    awk 'NR==FNR{keep[$1]=1;next} /^>/{id=$0; sub(/^>/,"",id); sub(/[[:space:]].*$/,"",id); print_flag=(id in keep)} { if(print_flag) print }' "$CID_TMP" "$FAS" > "$FAS_OUT"

    RPT="$ODIR/subset_report.${CAP}_${SAMPLE}.txt"
    TOTAL_ROWS=$(($(wc -l < "$CSV")-1))
    KEPT_ROWS=$(($(wc -l < "$CSV_OUT")-1))
    UBAR=$(awk -F, -v b="$BC_IDX" 'NR>1{t=$b; sub(/-1$/,"",t); print t}' "$CSV_OUT" | sort -u | wc -l | tr -d ' ')
    SEQS=$(grep -c '^>' "$FAS_OUT" || true)
    {
      echo "capture,$CAP"
      echo "sample,$SAMPLE"
      echo "input_csv_rows,$TOTAL_ROWS"
      echo "kept_csv_rows,$KEPT_ROWS"
      echo "unique_barcodes_kept,$UBAR"
      echo "fasta_contigs_kept,$SEQS"
    } > "$RPT"
    log "[VDJ] kept_csv_rows=$KEPT_ROWS contigs=$SEQS"
  done

  IND="$VDJROOT/indices/index_per_sample_vdj.csv"
  echo "capture,sample,out_dir,kept_csv_rows,fasta_contigs_kept,unique_barcodes_kept" > "$IND"
  find "$VDJROOT/per_sample" -type f -name "subset_report.*.txt" | sort | while read -r R; do
    CAP=$(awk -F, '$1=="capture"{print $2}' "$R")
    SMP=$(awk -F, '$1=="sample"{print $2}' "$R")
    ROWS=$(awk -F, '$1=="kept_csv_rows"{print $2}' "$R")
    FASTA=$(awk -F, '$1=="fasta_contigs_kept"{print $2}' "$R")
    UBAR=$(awk -F, '$1=="unique_barcodes_kept"{print $2}' "$R")
    ODIR="$(dirname "$R")"
    printf "%s,%s,%s,%s,%s,%s\n" "$CAP" "$SMP" "$ODIR" "$ROWS" "$FASTA" "$UBAR" >> "$IND"
  done
  log "VDJ summary: $IND"
else
  log "VDJ: No samples to subset."
fi

log "=== PIPELINE DONE ==="
log "Key outputs:"
log "• Seurat merged:        $MERGED_RDS"
log "• HTO CSV:              $HTO_CSV"
log "• HashSolo results:     $HASHROOT (incl. hashsolo_cells_full.csv)"
log "• Named Seurat:         $NAMED_RDS"
log "• QC counts (PNGs):     $QCCOUNTS/out"
log "• Confident singlets:   $CONF_RDS"
log "• Split per capture:    $SPLITROOT (indices/...)"
log "• VDJ per-sample out:   $VDJROOT (indices/index_per_sample_vdj.csv)"
