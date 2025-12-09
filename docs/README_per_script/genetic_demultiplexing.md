# Genetic Demultiplexing – GRO11679 (Souporcell + HashSolo)

This README describes how to use the two scripts for genetic demultiplexing with Souporcell, annotation with HashSolo, and plotting for project **GRO11679**.

Scripts:
- **Script 1 (Souporcell launcher):** `genetic_demultiplexing_script_1_GRO11679.sbatch`
- **Script 2 (Annotation + plotting):** `genetic_demultiplexing_script_2_annotation_and_plot__from_HS_GRO11679`

--------------------------------------------------------------------------------
## 1. Overview

1. **Run Souporcell on each capture**  
   Script 1 is a *launcher* run on the login node. It creates and submits a Slurm **array** job to run Souporcell for each capture (C1, C6, …).

2. **Annotate Souporcell clusters with HashSolo and plot**  
   Script 2 is submitted with `sbatch`. It:
   - Links Souporcell clusters to HashSolo singlet calls and patient IDs.
   - Writes annotated tables.
   - Generates per-capture barplots (cells per Souporcell cluster, coloured by patient).

--------------------------------------------------------------------------------
## 2. Requirements

- SLURM-based HPC cluster.
- Working **Souporcell** environment:
  - `MINIFORGE_ROOT` and `CONDA_ENV_PATH` with:
    - `souporcell_pipeline.py`
    - `samtools`
- Working **R** environment:
  - `R_BIN` pointing to an Rscript executable.
  - R packages: `readr`, `dplyr`, `stringr`, `purrr`, `ggplot2`.
- Completed upstream runs:
  - **Cell Ranger multi** runs per capture under `CR_RUN_ROOT`.
  - **HashSolo** run providing `hashsolo_cells_full.csv`.
  - A `captures.tsv` file mapping HTOs to sample IDs.

--------------------------------------------------------------------------------
## 3. Script 1 – Run Souporcell per capture

**File:** `genetic_demultiplexing_script_1_GRO11679.sbatch`  
**Run:** on the login node via `bash`, **not** with `sbatch`.

### 3.1. Key user settings to edit

Inside the script:

- Demultiplexing / project paths:
  - `DEMUX_ROOT="/hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing2"`
  - `PROJ="${DEMUX_ROOT}"` (output root for this workflow)
- Reference FASTA:
  - `REF_FASTA="/hpc/.../REF_GEX/refdata-gex-GRCh38-2024-A/fasta/genome.fa"`
- Cell Ranger multi root:
  - `CR_RUN_ROOT="/hpc/.../GRO11679/alignment_CR/cellranger_multi_2025-12-01"`
- Conda (Miniforge) + env with Souporcell:
  - `MINIFORGE_ROOT="/hpc/.../software/venvs4/miniforge3"`
  - `CONDA_ENV_PATH="/hpc/.../software/venvs4/souporcell"`
- Captures and K values:
  - `CAPTURE_IDS=(C1 C6)`     # captures to process
  - `K_VALUES=(2  2)`         # same length and order as CAPTURE_IDS
- Email notification (optional):
  - `MAIL_USER="${MAIL_USER:-egrossewichtrup@umcutrecht.nl}"`  
    (set empty to disable or override via env).

### 3.2. What Script 1 does

For all captures in `CAPTURE_IDS`:

- Creates project directories (if not present):
  - `${PROJ}`, `${PROJ}/demux`, `${PROJ}/tmp`, `${PROJ}/logs`, `${PROJ}/sbatch`
- Writes a timestamped Slurm **array script** to:
  - `${PROJ}/sbatch/run_souporcell_array_<timestamp>.sbatch`
- Submits an array job with one array task per capture.
- For each array task (capture `C`):
  - Finds the correct Cell Ranger multi directory for `C` under `CR_RUN_ROOT` using `find`.
  - Derives:
    - BAM:  
      `${MULTI_DIR}/outs/per_sample_outs/<multi_basename>/count/sample_alignments.bam`
    - Barcodes (gz):  
      `${MULTI_DIR}/outs/per_sample_outs/<multi_basename>/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz`
  - Creates:
    - Output root: `${PROJ}/demux/C`
    - Souporcell outdir: `${PROJ}/demux/C/C_souporcell`
    - Job-local temp dir under `${PROJ}/tmp`.
  - Activates the conda env (`CONDA_ENV_PATH`) and checks:
    - `souporcell_pipeline.py`
    - `samtools`
  - Indexes FASTA and BAM if needed (`samtools faidx`, `samtools index`).
  - Unzips barcodes to `${OUTDIR}/barcodes.tsv`.
  - Runs Souporcell:

    ```bash
    souporcell_pipeline.py \
      -i "${BAM}" \
      -b "${BARCODES_TSV}" \
      -f "${REF_FASTA}" \
      -t "${SLURM_CPUS_PER_TASK}" \
      -o "${OUTDIR}" \
      -k "${K}"
    ```

### 3.3. Inputs and outputs (Script 1)

**Inputs:**
- 10x reference FASTA: `REF_FASTA`.
- Cell Ranger multi outputs under `CR_RUN_ROOT` with:
  - Per-capture subdirs named like `C1 *` and `multi_C1_*`.
- Souporcell environment via `MINIFORGE_ROOT` and `CONDA_ENV_PATH`.

**Outputs (per capture `C`):**
- Souporcell output directory:
  - `${PROJ}/demux/C/C_souporcell/` (includes `clusters.tsv`).
- Slurm logs:
  - `${PROJ}/logs/souporcell_multi_<JOBID>_<TASKID>.out`
  - `${PROJ}/logs/souporcell_multi_<JOBID>_<TASKID>.err`

### 3.4. How to start Script 1

On the login node:

```bash
bash genetic_demultiplexing_script_1_GRO11679.sbatch
