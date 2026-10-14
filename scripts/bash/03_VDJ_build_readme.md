# JIA B-cell receptor (VDJ) processing — README v1.3

## 1. Input and output files / paths

**Project root (HPC):**

```text
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-10_VDJ_build_own_CR
```

**Inputs**
- **Sample sheet:** `00_config/samples_JIA_16.tsv` — lists the 16 samples, their patient/tissue identities, and the exact Cell Ranger input paths (`vdj_dir`).
- **Cell Ranger VDJ folders**, from **our own alignments**:
  - Blood (GRO11679): `/hpc/dla_lti/yermanos_group/GRO11679/2026-06_alignment_CR/cellranger_multi_outputs/`
  - Synovial fluid (GRO11735): `/hpc/dla_lti/yermanos_group/GRO11735/alignment_CR/`
- **Files inside each sample's `vdj_b` folder:** `filtered_contig_annotations.csv`, `filtered_contig.fasta`, `consensus_annotations.csv`, `consensus.fasta`, `clonotypes.csv`, and `concat_ref.fasta`.
- **Script:** `04_scripts/01_unified_jia_vdj_build_v1_3.sbatch` (generates a companion `.R` file in `04_scripts/`).

**Outputs** (all relative to the project root)

```text
02_results/vdj_build_v1_3/
├── VDJ_build_all_16_samples.rds.xz         ← MAIN DATASET: all cells together
├── by_sample/VDJ_build_<sample_id>.rds.xz  ← 16 individual sample tables
├── by_patient/VDJ_build_<patient_id>.rds.xz ← 13 patient-level tables
├── QC_counts_by_sample.tsv                 ← QC per sample
├── QC_counts_by_capture.tsv                ← QC per capture
├── sample_input_manifest.tsv               ← Input paths used
├── sample_mapping_used.tsv                 ← Saved copy of the sample sheet
├── sessionInfo.txt / R_source_used.R       ← Software versions and code
└── COMPLETION_REPORT.txt                   ← Created on successful completion
```

Job logs are written to `03_logs/JIA_VDJ_ownCR_<jobid>.out` and `.err`. A failed job keeps partial results and creates `FAILURE_REPORT.txt` (and possibly `diagnostics/ERROR_TRACE_<sample>.txt`).

## 2. What are we doing, and why?

**Cell Ranger has already aligned and identified the B-cell receptor (BCR) sequences.** This script does **not** redo the alignment. Instead, it uses the R package **Platypus** to organize those results into one consistent, cell-by-cell dataset. This gives us a starting point for studying B-cell clones and comparing repertoires across blood and joint fluid, and later linking BCR data to single-cell gene expression.

We process **16 samples from 13 people**: **11 blood** and **5 synovial-fluid** samples across **10 sequencing captures**.

## 3. What happens in the script?

1. **Check:** Verify sample identities, input files, and paths. Also make sure we are not accidentally using duplicated/pooled barcodes.
2. **Make formats compatible:** Apply temporary fixes so newer Cell Ranger outputs work with our installed Platypus/R packages. **The original Cell Ranger files remain unchanged.**
3. **Build BCR tables:** Run `VDJ_build()` for each sample. Its settings remove divergent cells, require complete BCR information, and enable germline trimming. Consequently, some input cells may be excluded.
4. **Label and combine:** Keep original Cell Ranger barcodes and clonotype IDs, add unambiguous sample/patient labels, save sample- and patient-level tables, and create the combined dataset and QC summaries.

**Important:** A clonotype label is kept specific to its original sample. The script does **not** prove that two cells from different samples or joints belong to the same biological clone. That needs a separate sequence-based comparison.

## 4. How to run / confirm success

```bash
ROOT="/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-10_VDJ_build_own_CR"
SCRIPT="$ROOT/04_scripts/01_unified_jia_vdj_build_v1_3.sbatch"
bash "$SCRIPT" --check    # Check inputs only; do not build results
sbatch "$SCRIPT"          # Submit full processing job
```

A successful run writes `02_results/vdj_build_v1_3/COMPLETION_REPORT.txt` with `STATUS: COMPLETE`, plus a `.COMPLETE` marker. The script refuses to overwrite an existing results folder: use a new version/output location for any rerun.
