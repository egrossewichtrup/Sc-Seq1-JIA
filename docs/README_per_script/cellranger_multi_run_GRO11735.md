TITLE: Cell Ranger multi — GRO11735 (GEX + VDJ-B, lanes 5–8)
VERSION: v0.1.0
LAST UPDATED: 2025-10-29
OWNER: Enno Große Wichtrup
CONTACT: e.a.grossewichtrup@umcutrecht.nl

PURPOSE
- Run cellranger multi for non-hashed human CD19+ B-cell samples (sorted from PBMCs of healthy controls and JIA patients) on GRCh38.
- Processes 5' Gene Expression (GEX) + BCR (VDJ-B) on lanes L005–L008.
- Builds the per-sample .multi.csv on the fly from fastq/pairs.tsv.

WHEN TO USE
- Batch processing of multiple GEX/VDJ sample pairs with a stable lane layout on Chromium X / GEM-X 5'.
- Dataset: human, B cells from PBMCs (healthy control + JIA), non-hashed, GRCh38.

SCRIPT AND PATHS
- Script file: cellranger_multi_run_GRO11735.sbatch
- Project root (default): /hpc/dla_lti/yermanos_group/GRO11735
- FASTQ directory: $ROOT/fastq
- Mapping file: $FASTQ_ROOT/pairs.tsv
- Outputs root (time-stamped): $ROOT/cellranger_<YYYY-MM-DD_HHMM>/

REFERENCES USED BY THE SCRIPT
- Gene Expression (GRCh38): /hpc/dla_lti/yermanos_group/REF_GEX/refdata-gex-GRCh38-2024-A
- VDJ-B (GRCh38): /hpc/dla_lti/yermanos_group/REF_VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

PLATFORM / RUN METADATA (INFO)
- Application: 10x Single Cell
- Platform: Chromium X
- Run Type: GEM-X 5' Gene Expression
- Add-ons: BCR (used), Cell surface protein (not used by this script)

----------------------------------------
QUICK START
----------------------------------------
1) Sanity checks
   ROOT="/hpc/dla_lti/yermanos_group/GRO11735"
   FASTQ_ROOT="$ROOT/fastq"
   PAIRS="$FASTQ_ROOT/pairs.tsv"
   ls -l "$PAIRS"; head -n 3 "$PAIRS"

2) Smoke test on the first GEX/VDJ pair (line 1)
   sbatch --array=1-1 cellranger_multi_run_GRO11735.sbatch

3) Full run for all pairs
   N=$(wc -l < "$PAIRS")
   sbatch --array=1-"$N" cellranger_multi_run_GRO11735.sbatch

4) Tail logs
   tail -f "$ROOT"/slurm/multi_11735_*.out "$ROOT"/slurm/multi_11735_*.err

----------------------------------------
WHAT THE SCRIPT DOES
----------------------------------------
- Reads line N from fastq/pairs.tsv (controlled by SLURM_ARRAY_TASK_ID).
  - Column 1 -> GEX_PREFIX
  - Column 2 -> VDJ_PREFIX
- Verifies FASTQs exist for lanes 5, 6, 7, 8 for both prefixes.
- Builds a per-sample .multi.csv with:
  * [gene-expression] reference = GRCh38 2024-A, create-bam = true
  * [vdj]             reference = GRCh38 alts ensembl 7.1.0 (VDJ-B)
  * [libraries]       two rows, feature_types = Gene Expression and VDJ-B, lanes = 5|6|7|8
- Runs cellranger multi in local mode using Slurm-provided cores and memory.
- Writes time-stamped results under "$ROOT/cellranger_<YYYY-MM-DD_HHMM>/<ID>/outs/".

ID FORMAT
- The script creates an ID like: GRO11735_<hash>_L5-8 (hash computed from GEX_PREFIX_L5-8, length-safe for Cell Ranger).

----------------------------------------
INPUTS (EXACT)
----------------------------------------
1) Mapping file (mandatory)
   Path: $FASTQ_ROOT/pairs.tsv
   Format: 2 columns, whitespace-separated (no header)
   Meaning: Column 1 = GEX prefix, Column 2 = VDJ-B prefix
   Example lines:
   P002-Sc-Seq1-LK-SFMC-Bcells-17-07-2025_A235V25LT3    P002-Sc-Seq1-LK-SFMC-Bcells-17-07-2025-BCR_A235V25LT3
   P002-Sc-Seq1-RK-SFMC-Bcells-17-07-2025_A235V25LT3    P002-Sc-Seq1-RK-SFMC-Bcells-17-07-2025-BCR_A235V25LT3
   P001-Sc-Seq1-LK-SFMC-Bcells-16-07-2025_A235V25LT3    P001-Sc-Seq1-LK-SFMC-Bcells-16-07-2025-BCR_A235V25LT3
   P004-Sc-Seq1-K-SFMC-Bcells-23-07-2025_A235V25LT3     P004-Sc-Seq1-K-SFMC-Bcells-23-07-2025-BCR_A235V25LT3
   P001-Sc-Seq1-RK-SFMC-Bcells-16-07-2025_A235V25LT3    P001-Sc-Seq1-RK-SFMC-Bcells-16-07-2025-BCR_A235V25LT3

2) FASTQ files (mandatory)
   Directory: $FASTQ_ROOT/
   Required pattern (both GEX and VDJ prefixes present):
   <PREFIX>_S<NN>_L00{5,6,7,8}_{I1|I2|R1|R2}_001.fastq.gz

   Example subset (abbreviated):
   P001-...-BCR_A235V25LT3_S1_L005_R1_001.fastq.gz
   P001-...-BCR_A235V25LT3_S1_L005_R2_001.fastq.gz
   P001-..._A235V25LT3_S2_L005_R1_001.fastq.gz
   P001-..._A235V25LT3_S2_L005_R2_001.fastq.gz
   ... (lanes L005–L008 for both GEX and BCR prefixes; include I1/I2/R1/R2)

3) References (mandatory; hard-coded in script)
   GEX: /hpc/dla_lti/yermanos_group/REF_GEX/refdata-gex-GRCh38-2024-A
   VDJ: /hpc/dla_lti/yermanos_group/REF_VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0

----------------------------------------
OUTPUTS (PER SAMPLE ID)
----------------------------------------
Root folder for a given run:
  $ROOT/cellranger_<YYYY-MM-DD_HHMM>/<ID>/outs/

Key outputs:
- Web summary (HTML)
  outs/per_sample_outs/<ID>/web_summary.html

- Filtered GEX matrix (HDF5)
  outs/per_sample_outs/<ID>/count/sample_filtered_feature_bc_matrix.h5

- BCR contig annotations (CSV)
  outs/per_sample_outs/<ID>/vdj_b/filtered_contig_annotations.csv

- Alignments (BAM and index)
  outs/per_sample_outs/<ID>/count/sample_alignments.bam
  outs/per_sample_outs/<ID>/count/sample_alignments.bam.bai

Other useful artifacts:
- outs/metrics_summary.csv (pipeline-level metrics)
- <ID>/_log (Martian pipeline logs and timeline)

----------------------------------------
SLURM PROFILE (AS IN SCRIPT)
----------------------------------------
--partition=cpu
--account=dla_lti
--cpus-per-task=32
--mem=96G
--time=24:00:00
--output=/hpc/dla_lti/yermanos_group/GRO11735/slurm/%x_%A_%a.out
--error=/hpc/dla_lti/yermanos_group/GRO11735/slurm/%x_%A_%a.err

CELL RANGER RUNTIME PARAMS
- --jobmode=local
- --localcores=$SLURM_CPUS_PER_TASK (32 by default)
- --localmem=90 (GB)
- --disable-ui
- create-bam=true

----------------------------------------
CONFIGURATION KNOBS (EDIT AT TOP OF SCRIPT)
----------------------------------------
CR             = /hpc/dla_lti/yermanos_group/software/cellranger-9.0.1/bin/cellranger
ROOT           = /hpc/dla_lti/yermanos_group/GRO11735
FASTQ_ROOT     = $ROOT/fastq
PAIRS          = $FASTQ_ROOT/pairs.tsv
lanes          = 5|6|7|8
VDJ feature    = VDJ-B
create-bam     = true

Note on chemistry:
- Autodetect should be fine for Chromium X GEM-X 5' GEX + BCR.
- If you must override chemistry, add the appropriate flag or CSV field as per Cell Ranger docs.

----------------------------------------
CHECKS AND VALIDATION
----------------------------------------
- FASTQ presence check: the script confirms at least one FASTQ exists for each of lanes 5–8 for both prefixes.
- CSV hygiene: removes any accidental "chain" column; forces VDJ-B.
- Stale artifacts: removes any prior same-ID folders within the current run root.

Quick validation workflow:
1) Submit a one-line array:
   sbatch --array=1-1 cellranger_multi_run_GRO11735.sbatch
2) Inspect the constructed CSV echoed in the .out:
   grep -n "=== CSV FINAL ===" "$ROOT"/slurm/multi_11735_*_*.out | tail -n1 | cut -d: -f1 | xargs -I{} sed -n '1,200p' {}

----------------------------------------
TROUBLESHOOTING
----------------------------------------
Symptom: "ERROR: No line found in pairs.tsv"
Cause: Array index exceeds number of lines
Fix: Use N=$(wc -l < pairs.tsv) and submit with --array=1-$N

Symptom: FASTQ check fails for a prefix/lane
Cause: Wrong prefixes in pairs.tsv or missing L005–L008 files
Fix: Verify pairs.tsv lines and confirm files exist for each lane

Symptom: Slow run or out-of-memory
Cause: Under-provisioned resources
Fix: Increase --mem, --cpus-per-task, and the cellranger --localmem value

Symptom: Web summary missing
Cause: Early pipeline failure
Fix: Check slurm/*.err and <ID>/_log under outs; re-run after fixing root cause

----------------------------------------
REPRODUCIBILITY
----------------------------------------
- Record commit: git rev-parse HEAD >> run.log
- Record Cell Ranger version: "$CR" --version | sed 's/^/cellranger_version=/' >> run.log
- Record command: echo "cmd=cellranger multi --id=<ID> --csv=<CSV> ..." >> run.log
- Optionally export environment: mamba env export -n <env> > env.yaml  (or pip freeze > requirements.txt)

----------------------------------------
FILE MAP
----------------------------------------
./
├─ cellranger_multi_run_GRO11735.sbatch
├─ README.md
└─ fastq/
   ├─ pairs.tsv
   └─ <FASTQs for GEX + VDJ-B across L005–L008>

----------------------------------------
EXTENSIONS
----------------------------------------
- To include Cell Surface Protein (Feature Barcoding), add a third [libraries] row in the .multi.csv with feature_types=Antibody Capture and the appropriate CSP fastq_id pointing to the CSP prefix. Adjust references or feature refs as needed.
- For TCR (VDJ-T or VDJ-GD), change the VDJ reference block and the feature_types accordingly in the [libraries] section.

END OF FILE
