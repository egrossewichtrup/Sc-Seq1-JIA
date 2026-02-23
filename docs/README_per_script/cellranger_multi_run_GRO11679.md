TITLE: Cell Ranger multi — GRO11679 (GEX + CSP/HTO + VDJ-B, two runs A235V25LT3 & B22WFC2LT3, array sbatch)
VERSION: v0.2.0
LAST UPDATED: 2025-10-30
OWNER: Enno Große Wichtrup
CONTACT: e.a.grossewichtrup@umcutrecht.nl

PURPOSE
- Run Cell Ranger “multi” for the GRO11679 dataset (human, B cells from PBMCs; healthy control + JIA), combining:
  • 5′ Gene Expression (GEX)  • Cell Surface Protein / HTO (Antibody Capture)  • VDJ-B (BCR)
- Inputs are split across two sequencing runs (A235V25LT3 and B22WFC2LT3).
- One Slurm ARRAY job launches one task per capture; each task builds its own multi-config and runs `cellranger multi`.

WHEN TO USE
- You have multiple captures with the same library types (GEX + CSP/HTO + VDJ-B) across two runs and you want consistent per-capture outputs.
- Dataset specifics: human GRCh38; non-hashed for GEX/VDJ (HTOs used for CSP only).

SCRIPT & PATHS
- Script file name: cellranger_multi_run_GRO11679.sbatch
- Recommended location: /hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/scripts/cellranger_multi_run_GRO11679.sbatch
- Paths inside the script are absolute; you can submit from anywhere with `sbatch /hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/scripts/cellranger_multi_run_GRO11679.sbatch`.

PLATFORM / RUN META
- Application: 10x Single Cell  • Platform: Chromium X  • Run Type: GEM-X 5′ GEX
- Libraries: GEX (Gene Expression), Antibody Capture (CSP / HTO), VDJ-B (BCR)
- Cell Ranger: v9.0.1  • Binary used in script: /hpc/dla_lti/yermanos_group/software/cellranger-9.0.1/bin

REFERENCES (as used by the script)
- GEX (GRCh38): /hpc/dla_lti/yermanos_group/REF_GEX/refdata-gex-GRCh38-2024-A
- VDJ-B (GRCh38): /hpc/dla_lti/yermanos_group/REF_VDJ/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0
- Feature barcode CSV: $PROJ/feature_ref/featurebarcode.csv  (created below; required for HTOs)

----------------------------------------
QUICK START (COPY–PASTE)
----------------------------------------
0) Set your analysis root and make the feature-ref folder
   export PROJ="/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348"
   mkdir -p "$PROJ/feature_ref"
   # Note: the script itself will also create "$PROJ/logs" and a dated run root.

1) Create / verify the feature barcode CSV (WITH HEADER, HTO1–HTO4)
cat > "/hpc/dla_lti/yermanos_group/GRO11679/alignment_CR/feature_ref/featurebarcode.csv" <<'CSV'
id,name,read,pattern,sequence,feature_type
HTO1,HTO1,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
HTO2,HTO2,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
HTO3,HTO3,R2,5PNNNNNNNNNN(BC),TTCCGCCTCTCTTTG,Antibody Capture
HTO4,HTO4,R2,5PNNNNNNNNNN(BC),AGTAAGTTCAGCGTA,Antibody Capture
CSV

2) Save the array script
   nano "$PROJ/scripts/cellranger_multi_run_GRO11679.sbatch"
   # Use the sbatch header exactly as below (matches the current script):
   #SBATCH --job-name=cr_multi_array
   #SBATCH --time=2-00:00:00
   #SBATCH --cpus-per-task=16
   #SBATCH --mem=128G
   #SBATCH --output=%x_%A_%a.out
   #SBATCH --error=%x_%A_%a.err
   #SBATCH --array=1-5
   # (Optional) If you prefer logs inside $PROJ/logs, change the two lines to:
   #SBATCH --output=/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/logs/%x_%A_%a.out
   #SBATCH --error=/hpc/dla_lti/yermanos_group/GRO11679/analysis_20250925_140348/logs/%x_%A_%a.err

3) Submit the array (5 captures)
   sbatch "$PROJ/scripts/cellranger_multi_run_GRO11679.sbatch"

4) Monitor
   squeue -u "$USER"
   # By default (as in the script), Slurm logs are written to the directory you ran sbatch from:
   tail -f cr_multi_array_*_*.out cr_multi_array_*_*.err
   # If you changed the SBATCH lines to use $PROJ/logs:
   tail -f "$PROJ"/logs/cr_multi_array_*_*.out "$PROJ"/logs/cr_multi_array_*_*.err

----------------------------------------
WHAT THE SCRIPT DOES (ARRAY MODE)
----------------------------------------
- Defines a CAPTURE list (exact labels used for output folder names):
  C1 ScSeq1-Capture1-HC001-HC002-10-07-2025
  C6 ScSeq1-Capture6-P003-HC003-21-07-2025
  C8 ScSeq1-Capture8-P005-P006-24-07-2025
  C9 ScSeq1-Capture9-P002blood-P007-P008-30-07-2025
  C10 ScSeq1-Capture10-P009-P010-31-07-2025
- Creates a dated run root on the fly: $PROJ/cellranger_multi_<YYYY-MM-DD>
- For the current array task, builds a per-capture CSV with:
  [gene-expression] (GRCh38, create-bam=false)
  [vdj] (GRCh38 VDJ-B ref)
  [feature] (HTO feature_ref CSV)
  [libraries] rows that reference BOTH runs for each library type:
    <stub>_A235V25LT3, <FASTQ_A>, Gene Expression
    <stub>_B22WFC2LT3, <FASTQ_B>, Gene Expression
    <stub>-CSP_A235V25LT3, <FASTQ_A>, Antibody Capture
    <stub>-CSP_B22WFC2LT3, <FASTQ_B>, Antibody Capture
    <stub>-BCR_A235V25LT3, <FASTQ_A>, VDJ-B
    <stub>-BCR_B22WFC2LT3, <FASTQ_B>, VDJ-B
- Runs: `cellranger multi --id=multi_<safe_label> --csv=<cfg>`
  (safe_label = label with spaces replaced by underscores)

----------------------------------------
INPUTS (EXACT)
----------------------------------------
A) FASTQ directories (two runs):
   Run A: /hpc/dla_lti/yermanos_group/GRO11679/fastq/A235V25LT3
   Run B: /hpc/dla_lti/yermanos_group/GRO11679/fastq/B22WFC2LT3

B) FASTQ name patterns (per library type; omit S## in fastq_id)
   <stub>_A235V25LT3_SNN_L00{5..8}_{I1|I2|R1|R2}_001.fastq.gz
   <stub>_B22WFC2LT3_SNN_L00{5..8}_{I1|I2|R1|R2}_001.fastq.gz
   <stub>-CSP_A235V25LT3_SNN_L00{5..8}_{I1|I2|R1|R2}_001.fastq.gz
   <stub>-CSP_B22WFC2LT3_SNN_L00{5..8}_{I1|I2|R1|R2}_001.fastq.gz
   <stub>-BCR_A235V25LT3_SNN_L00{5..8}_{I1|I2|R1|R2}_001.fastq.gz
   <stub>-BCR_B22WFC2LT3_SNN_L00{5..8}_{I1|I2|R1|R2}_001.fastq.gz

   Examples (truncated):
   P001-...-BCR_A235V25LT3_S1_L005_R1_001.fastq.gz
   P001-..._A235V25LT3_S2_L005_R2_001.fastq.gz
   … across L005–L008 for both runs and all three library types.

C) Feature barcode CSV (required for HTOs)
   Path: $PROJ/feature_ref/featurebarcode.csv
   Header: id,name,read,pattern,sequence,feature_type
   Entries: HTO1–HTO4 as listed in Quick Start step 1 (R2; pattern 5PNNNNNNNNNN(BC))

D) Script file
   $PROJ/scripts/cellranger_multi_run_GRO11679.sbatch

----------------------------------------
OUTPUTS (PER CAPTURE)
----------------------------------------
Run root (auto-dated):
  $PROJ/cellranger_multi_<YYYY-MM-DD>/

Per-capture working directory (keeps your full label with spaces):
  $PROJ/cellranger_multi_<YYYY-MM-DD>/<CAPTURE FULL LABEL>/

Cell Ranger outputs:
  $PROJ/cellranger_multi_<YYYY-MM-DD>/<CAPTURE FULL LABEL>/outs/

Key files inside outs/ (paths shown with <id>=multi_<safe_label>):
- Web summary (QC): outs/per_sample_outs/<id>/web_summary.html
- GEX filtered matrix (HDF5): outs/per_sample_outs/<id>/count/sample_filtered_feature_bc_matrix.h5
- VDJ-B contigs (CSV): outs/per_sample_outs/<id>/vdj_b/filtered_contig_annotations.csv
- Metrics: outs/metrics_summary.csv (+ per-library JSONs)
- (No BAMs; create-bam=false to save space)

Config used:
  $PROJ/cellranger_multi_<YYYY-MM-DD>/multi_config_<safe_label>.csv

Logs:
- Slurm (default with current SBATCH): %x_%A_%a.out and %x_%A_%a.err in the directory you ran `sbatch` from.
- (Optional) If you edited SBATCH to point to $PROJ/logs, they’ll be there.
- Cell Ranger runtime: outs/cellranger.log

----------------------------------------
FASTQ_ID NAMING RULES (CRITICAL)
----------------------------------------
✅ Use the prefix BEFORE `_SNN` and include the library tag exactly:
   <stub>_A235V25LT3
   <stub>_B22WFC2LT3
   <stub>-CSP_A235V25LT3
   <stub>-CSP_B22WFC2LT3
   <stub>-BCR_A235V25LT3
   <stub>-BCR_B22WFC2LT3

❌ Do NOT include `_SNN` in fastq_id values (Cell Ranger matches stems; `_SNN` is variable).
❌ Do NOT drop the `-CSP` or `-BCR` markers for those libraries.

----------------------------------------
SLURM RESOURCES (IN SCRIPT)
----------------------------------------
#SBATCH --job-name=cr_multi_array
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-5

Cell Ranger local mode settings are derived from SBATCH (localcores=cpus-per-task, localmem≈mem).

----------------------------------------
MONITORING & VERIFICATION
----------------------------------------
- Queue:       squeue -u "$USER"
- Slurm logs (default):  tail -f cr_multi_array_*_*.out
- CR log (example capture):  tail -n 120 "$PROJ"/cellranger_multi_*/"C1 ScSeq1-Capture1-HC001-HC002-10-07-2025"/outs/cellranger.log

Quick sanity (config for one capture):
- Pick a config: ls "$PROJ"/cellranger_multi_*/multi_config_C1_ScSeq1-Capture1-*.csv | tail -n1 | xargs -I{} sed -n '1,120p' {}

----------------------------------------
COMMON ERRORS & FIXES
----------------------------------------
“Feature reference file header does not contain required fields”
→ Ensure the header line exactly matches:
   id,name,read,pattern,sequence,feature_type

“Requested sample(s) not found in fastq directory … Available samples: …”
→ fastq_id mismatch. Remove `_SNN`; keep `-CSP` / `-BCR` where appropriate.
→ Compare your config fastq_id values to the stems printed in the error.

Preflight passes but pipeline fails early
→ Check outs/cellranger.log and the per-library JSONs under outs.

----------------------------------------
RERUNS / CLEANUP
----------------------------------------
- Cell Ranger won’t overwrite an existing `--id`. Either bump the id or remove the prior pipestance folder.
- Per-capture cleanup example (uses safe labels):
  for L in "C1 ScSeq1-Capture1-HC001-HC002-10-07-2025" \
           "C6 ScSeq1-Capture6-P003-HC003-21-07-2025" \
           "C8 ScSeq1-Capture8-P005-P006-24-07-2025" \
           "C9 ScSeq1-Capture9-P002blood-P007-P008-30-07-2025" \
           "C10 ScSeq1-Capture10-P009-P010-31-07-2025"; do
    SAFE="${L// /_}"
    rm -rf "$PROJ"/cellranger_multi_*/"$L"           # work dir (quoted, has spaces)
    rm -f  "$PROJ"/cellranger_multi_*/multi_config_"$SAFE".csv
  done

----------------------------------------
CUSTOMIZATION
----------------------------------------
- Add/remove captures: edit the CAPS array in the sbatch and adjust `--array=1-N`.
- Resources: tune `--cpus-per-task` / `--mem` to cluster norms.
- Chemistry: autodetect is fine for GEM-X 5′; override only if needed.

----------------------------------------
FILE MAP (TYPICAL AFTER A RUN)
----------------------------------------
$PROJ/
├─ scripts/
│  └─ cellranger_multi_run_GRO11679.sbatch
├─ feature_ref/
│  └─ featurebarcode.csv
├─ cellranger_multi_<YYYY-MM-DD>/
│  ├─ multi_config_C1_ScSeq1-Capture1-HC001-HC002-10-07-2025.csv
│  ├─ multi_config_C6_ScSeq1-Capture6-P003-HC003-21-07-2025.csv
│  ├─ …
│  ├─ C1 ScSeq1-Capture1-HC001-HC002-10-07-2025/
│  │  └─ outs/ (web_summary.html, per_sample_outs/<id>/{count,vdj_b}, metrics files, logs)
│  ├─ C6 ScSeq1-Capture6-P003-HC003-21-07-2025/
│  │  └─ outs/ …
│  └─ …
# Note: With the current SBATCH lines, Slurm logs (%x_%A_%a.out/.err) are in the directory you ran `sbatch` from.
#       If you edit SBATCH to use $PROJ/logs, they will instead appear under $PROJ/logs/.

END OF FILE
