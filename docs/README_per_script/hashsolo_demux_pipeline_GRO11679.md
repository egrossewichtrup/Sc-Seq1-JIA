DEMUX_PIPELINE — Hashing-based demultiplexing and VDJ subsetting on SLURM
=========================================================================

Overview
--------
This SLURM-driven bash pipeline merges multiple 10x Cell Ranger multi runs, performs HashSolo demultiplexing from HTO counts, attaches demux metadata back to Seurat, generates QC summaries, splits confident singlets per capture and per sample, and subsets VDJ contigs per sample. It is designed for HPC clusters with SLURM and pre-existing Miniforge/R environments.

Main Features
-------------
1) Robust discovery of Cell Ranger outputs even when folders contain spaces.  
2) Seurat merge of Gene Expression matrices and export of raw Antibody Capture (HTO) counts.  
3) HashSolo demultiplexing per capture with QC PDF reports and a combined CSV of cell-level calls.  
4) Seurat object augmented with HS_* metadata and mapped to short/full sample names from a user config.  
5) PNG QC plots of confident singlet counts and fractions per capture.  
6) Split of confident singlets per capture and per sample, with indices and barcode lists.  
7) VDJ subsetting per sample to produce filtered CSV and FASTA plus a manifest and summary index.  

Prerequisites
-------------
**Cluster:** SLURM scheduler  
**File system:** Access to your Cell Ranger multi outputs (`outs/per_sample_outs/...`)  
**Environments:**  
- Miniforge conda at `VENV_ROOT/venvs/miniforge3` with `conda.sh`  
- Python env at `VENV_ROOT/venvs/py` containing: `anndata`, `scanpy`, `scanpy-external`, `numpy`, `pandas`, `matplotlib`  
- R at `VENV_ROOT/venvs/r/bin/Rscript` with packages: `Seurat`, `Matrix`, `hdf5r`, `readr`, `dplyr`, `tibble`, `stringr`, `ggplot2`, `tidyr`  

Inputs and Assumptions
----------------------
- `CR_ROOT` points to a dated Cell Ranger multi root containing per-capture folders named like `multi_C1_*`, each with `outs/per_sample_outs/<runfolder>/count/sample_filtered_feature_bc_matrix.h5` and optionally `outs/per_sample_outs/<runfolder>/vdj_b/{filtered_contig_annotations.csv, filtered_contig.fasta}`.  
- The pipeline creates a self-contained run directory under `PROJ/demux` with time-stamped subfolders.  
- You provide a per-capture TSV `captures.tsv` (the script writes your initial table and moves it into the run folder automatically).  

Quick Start
-----------
1. **Edit the USER SETTINGS block** at the top of the script:
   - `PROJ`: project root for this run (outputs will be created beneath this)  
   - `VENV_ROOT`: root that contains `venvs/miniforge3`, `venvs/py`, and `venvs/r`  
   - `CR_ROOT`: root of your Cell Ranger multi runs for this batch  
   - Optional: adjust `HASH_SOLO_PRIORS` and `SINGLET_CONFIDENCE_MIN`

2. **Edit the CAPTURE CONFIG** in the script to match your data.  
   One line per capture with tab-separated columns:  
   `cap	runfolder	htos	map_short	map_full`

   Example lines:
   ```
   C1	multi_C1_ScSeq1-Capture1-HC001-HC002-10-07-2025	HTO1,HTO2	HTO1=HC001;HTO2=HC002	HTO1=HC001_Sc-Seq1_PBMC_Bcells_10.07.2025;HTO2=HC002_Sc-Seq1_PBMC_Bcells_10.07.2025
   C6	multi_C6_ScSeq1-Capture6-P003-HC003-21-07-2025	HTO1,HTO2	HTO1=P003;HTO2=HC003	HTO1=P003_Sc-Seq1_PBMC_Bcells_21.07.2025;HTO2=HC003_Sc-Seq1_PBMC_Bcells_21.07.2025
   ```
   *Note:* HTO feature names must match exactly those in the Antibody Capture matrix.

3. **Submit to SLURM:**
   ```
   sbatch demux_pipeline.sh
   ```

SLURM Defaults
--------------
```
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your@email
#SBATCH --chdir=/hpc/dla_lti/yermanos_group/GRO11679/demu_script_test
```
Output and error files: `DEMUX_PIPELINE_<jobid>.out` and `.err`

What the Pipeline Does
----------------------
**Step 0: Preflight**  
Checks Rscript, conda.sh, and Python envs.

**Step 1 — Seurat merge (R)**  
Merges H5 matrices into a combined Seurat object and exports raw HTO counts.  
Outputs under `SEUROOT`:
- `combined_multi_seurat_clean_<timestamp>_<caps>.rds`  
- `hto_counts_raw_<timestamp>_<caps>.csv`  
- `logs/qc_summary_<tag>.txt`

**Step 2 — HashSolo demux (Python)**  
Runs HashSolo per capture, saving QC PDFs and combined CSV.  
Outputs under `HASHROOT`:
- `hto_counts_raw.h5ad`  
- `<cap>/hashsoloed_<cap>.h5ad`  
- `<cap>/hashsolo_qc_<cap>.pdf`  
- `hashsolo_cells_full.csv`

**Step 3 — Attach & map (R)**  
Adds HS_* metadata to Seurat and maps HTO calls to sample names.  
Outputs:
- `<...>_with_hashsolo.rds`  
- `<...>_with_hashsolo_named.rds`

**Step 4 — QC plots (R)**  
Creates stacked count and fraction plots.  
Outputs in `QCCOUNTS/out`:
- `stacked_counts_per_capture.png`  
- `stacked_fractions_per_capture.png`  
- `counts_per_capture.csv`  
- `counts_per_capture_sample.csv`

**Step 5 — Filter confident singlets (R)**  
Removes negatives/doublets, keeping confident singlets only.  
Outputs in `FILTERROOT`:
- `seurat_confident_singlets.rds`  
- `counts_by_capture_flag.csv`

**Step 6 — Split per capture/sample (R)**  
Creates per-capture and per-sample Seurat RDS, barcodes, and indices.  
Outputs in `SPLITROOT`:
- `seurat/` and `barcodes/` folders  
- `indices/index_per_capture.csv`  
- `indices/index_per_sample.csv`  
- `indices/counts_per_capture.csv`  
- `indices/counts_per_capture_sample.csv`

**Step 7 — VDJ subsetting per sample (bash)**  
Subsets VDJ CSV and FASTA per sample using barcode lists.  
Outputs in `VDJROOT`:
- `per_sample/<cap>_<sample>/filtered_contig_annotations.<cap>_<sample>.csv`  
- `per_sample/<cap>_<sample>/filtered_contig.<cap>_<sample>.fasta`  
- `per_sample/<cap>_<sample>/subset_report.<cap>_<sample>.txt`  
- `indices/step7_manifest.csv`  
- `indices/index_per_sample_vdj.csv`

Configuration Knobs
-------------------
| Variable | Description | Default |
|-----------|--------------|----------|
| `PROJ` | Base output directory | user-defined |
| `VENV_ROOT` | Root with envs (miniforge3, py, r) | user-defined |
| `CR_ROOT` | Cell Ranger multi root | user-defined |
| `R_BIN` | Path to Rscript | `$VENV_ROOT/venvs/r/bin/Rscript` |
| `HASH_SOLO_PRIORS` | Priors for HashSolo | `"0.01,0.80,0.19"` |
| `SINGLET_CONFIDENCE_MIN` | Min. singlet prob | `0.90` |

Key Output Locations (relative to PROJ)
---------------------------------------
- `demux/<timestamp>_demux_<caps>/` → run root, logs, and generated scripts  
- `seurat/<timestamp>_seurat_merge_<caps>/` → merged Seurat and HTO CSV  
- `hashsolo/<timestamp>_hashsolo_<caps>/` → HashSolo results  
- `qc/step4_1_counts_<timestamp>/out/` → QC PNGs and counts  
- `qc/step4_2_filter_<timestamp>/` → confident singlets  
- `seurat/<timestamp>_step4_3_split/` → per-capture/sample splits  
- `vdj/<timestamp>_VDJ_PER_SAMPLE/` → per-sample VDJ subsets and indices  

Troubleshooting
---------------
**ERROR: No Rscript found**  
→ Ensure `R_BIN` points to a valid Rscript.

**ERROR: Miniforge conda.sh not found**  
→ Check `VENV_ROOT/venvs/miniforge3/etc/profile.d/conda.sh`.

**ERROR: Python env not found**  
→ Verify `VENV_ROOT/venvs/py` exists with dependencies.

**FATAL: merged RDS or HTO CSV missing**  
→ Check that capture names match actual folder names in `outs/per_sample_outs`.  
→ Ensure Antibody Capture data exists.

**Missing or incomplete VDJ**  
→ The pipeline skips missing ones automatically.

**Missing HTO features**  
→ Confirm `htos` column matches Antibody Capture feature names exactly.

Reproducibility and Re-runs
---------------------------
- Each run creates a new time-stamped folder.  
- Original inputs remain unchanged.  
- To re-run: edit settings and resubmit; outputs will land in a new folder.

Citations / Attributions
------------------------
- Seurat (Satija Lab)  
- HashSolo (Scanpy external)  
- 10x Genomics Cell Ranger  

License
-------
Add your project’s license here.
