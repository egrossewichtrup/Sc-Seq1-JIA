# Genetic demultiplexing and patient annotation (GRO11679)

Here I described three SLURM-ready bash scripts to run a full genetic demultiplexing workflow for 10x single-cell data using **Souporcell**, link genetic clusters to **hash-tag oligo (HTO)** patient IDs, and generate QC barplots.

1. genetic_demultiplexing_script_1_GRO11679.sbatch
2. genetic_demultiplexing_script_2_patient_annotation_GRO11679.sbatch
3. plot_souporcell_barplot.sbatch

The scripts are written for a GRO11679-style project structure on an HPC, but paths and parameters can be adapted to other setups.

---

## Overview of scripts

### 1. `genetic_demultiplexing_script_1_GRO11679.sbatch`

Launcher script that:

- Locates the 10x BAM + barcodes for a given capture (`CAPTURE_ID`)
- Submits a **Souporcell** job via an auto-generated `sbatch` file
- Organises all outputs under:

    ${PROJ}/demux/genetic_demultiplexing_<timestamp>/
        └─ <CAPTURE_ID>_souporcell/

Key settings inside the script:

- DEMUX_ROOT / PROJ: root directory for this genetic demultiplexing project  
- REF_FASTA: GRCh38 reference FASTA for Souporcell  
- CR_RUN_ROOT: Cell Ranger multi BAM root (multi_<CAPTURE_ID>_bam)  
- MINIFORGE_ROOT / CONDA_ENV_PATH: conda/Miniforge env with Souporcell + samtools  
- CAPTURE_ID: capture to process (default C1)  
- K: number of donors (clusters) for Souporcell (default 2)  

How to run (login node):

    # Default: run for C1 with K=2
    bash genetic_demultiplexing_script_1_GRO11679.sbatch

    # Example: run for C2 with K=3
    CAPTURE_ID=C2 K=3 bash genetic_demultiplexing_script_1_GRO11679.sbatch

---

### 2. `genetic_demultiplexing_script_2_patient_annotation_GRO11679.sbatch`

SLURM job that:

- Auto-detects the latest HashSolo output:
  - `hashsolo_cells_full.csv` under `${DEMUX_ROOT}/hashsolo/...`
  - Corresponding `captures.tsv` (HTO → sample mapping)
- Auto-detects the latest Souporcell run for C1 under `${GENETIC_ROOT}`
- Joins both sources:
  - Uses HashSolo singlets and the HTO sample mapping to infer patient IDs per genetic cluster
  - Matches HashSolo barcodes to Souporcell barcodes via the 10x barcode core
- Writes per-capture outputs:

    <SP_ROOT>/souporcell_with_patient_<timestamp>/<CAPTURE_ID>/
        ├─ <CAPTURE_ID>_cluster_to_patient.tsv
        └─ <CAPTURE_ID>_clusters_with_patient.tsv

Where:

- `*_cluster_to_patient.tsv` maps Souporcell cluster IDs (`cluster_id`) to `patient_id_from_HTO`
- `*_clusters_with_patient.tsv` is the original `clusters.tsv` plus a `patient_id_from_HTO` column per barcode

Important user-configurable variables at the top:

- PROJ_MAIN, DEMUX_ROOT, GENETIC_ROOT  
- HS_CSV and SP_ROOT (optional overrides; otherwise auto-detected)  
- CAPTURES (currently `"C1"`)  
- CONF_T (HashSolo singlet probability threshold, default `0.90`)  
- R_BIN (Rscript binary, e.g. an R venv)  

How to run:

    sbatch genetic_demultiplexing_script_2_patient_annotation_GRO11679.sbatch

After completion, check:

    ls /hpc/dla_lti/yermanos_group/GRO11679/genetic_demultiplexing/demux/
    # browse to: genetic_demultiplexing_<DATE>_.../souporcell_with_patient_<DATE>/C1/

---

### 3. `plot_souporcell_barplot.sbatch`

SLURM job that:

- Reads the annotated clusters file:  
  `<SP_ROOT>/<ANNOT_SUBDIR>/<CAPTURE>/_clusters_with_patient.tsv`
- Filters to Souporcell singlets with a non-NA `patient_id_from_HTO`
- Produces two barplots per capture:
  - Absolute cell counts per Souporcell cluster, stacked by patient
  - Fractions per cluster (relative composition), stacked by patient

Output files (for `CAPTURE="C1"`):

- `C1_souporcell_barplot_cells_per_cluster_by_patient.png`
- `C1_souporcell_barplot_fraction_per_cluster_by_patient.png`

Configurable variables:

- CAPTURE (e.g. `"C1"`)  
- SP_ROOT (Souporcell demux root)  
- ANNOT_SUBDIR (e.g. `souporcell_with_patient_20251122_095429`)  
- R_BIN (Rscript binary)  

How to run:

    sbatch plot_souporcell_barplot.sbatch

Resulting PNGs are written next to the `*_clusters_with_patient.tsv` file for that capture.

---

## Dependencies

- SLURM
- 10x Cell Ranger multi outputs with:
  - BAM + `barcodes.tsv.gz` for each capture
- Souporcell installed in the conda env specified by `CONDA_ENV_PATH`
- `samtools` available in the same env
- R with the packages:
  - `readr`, `dplyr`, `stringr`, `purrr`, `ggplot2`

Paths and environment variables in the scripts should be adapted to your own HPC and project layout before use.
