README
VDJ_build on capture level, then split into per sample and per patient objects

Project root
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline


Goal
Create Platypus VDJ dataframes with the recommended VDJ_build settings
remove.divergent.cells TRUE
complete.cells.only TRUE
trim.germlines TRUE

Because some samples were demultiplexed and lacked consensus files, we ran VDJ_build on the full capture level Cell Ranger VDJ output directories, then split the resulting capture level VDJ tables into per sample objects using demultiplexed barcode lists. Finally we merged samples per patient into one patient level VDJ table.


Inputs

1) Capture level full Cell Ranger VDJ folders
Configured in
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/00_config/captures_full_vdj_dirs.tsv

These folders must contain at least
filtered_contig_annotations.csv
filtered_contig.fasta
consensus_annotations.csv
consensus.fasta
clonotypes.csv
concat_ref.fasta

In our setup the capture level folders for C1, C6, C8, C9, C10 are the cellranger_multi per sample outs vdj_b folders under
/hpc/dla_lti/yermanos_group/GRO11679/alignment_CR/cellranger_multi_2025-12-01/...

Other captures can point to core data VDJ folders when those contain full Cell Ranger outputs.


2) Demultiplexed barcode lists per sample
Configured in
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/00_config/samples_demux_barcodes.tsv

For each sample the script uses
demux_filtered_contig_annotations_csv

Example path pattern
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2025-12_core_data/<capture>/VDJ/<capture>_<sample>/filtered_contig_annotations.csv

The barcode column in these files defines which 10x cell barcodes belong to the sample. Those barcodes are used to split the capture level VDJ table into per sample VDJ tables.


3) Reference core data used for demux barcode files
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2025-12_core_data/


How it was run

Slurm script
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/scripts/01_vdj_build_capture_then_split.sbatch

It creates and runs an R script at
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/04_scripts/01_vdj_build_capture_then_split.R

Environment
Conda R environment
/hpc/dla_lti/yermanos_group/software/venvs4/r


What the pipeline does

Step A
Run Platypus VDJ_build once per capture using the capture full VDJ directory.

Step B
Prefix clonotype_id with capture id to avoid clashes across captures.
A copy of the original clonotype id is kept in clonotype_id_10x_raw.

Step C
Split each capture level VDJ table into per sample VDJ tables by matching cell barcodes to the demux filtered_contig_annotations.csv barcode list.

Step D
Merge per sample VDJ tables into one per patient VDJ table by row binding all samples from that patient.


Outputs

All outputs live under
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/02_results/vdj_build/

1) Capture level VDJ tables
by_capture/<capture_id>/VDJ_build_capture_<capture_id>.rds.xz

Example
.../by_capture/C1/VDJ_build_capture_C1.rds.xz

2) Sample level VDJ tables
by_sample/<sample_id>/VDJ_build_<sample_id>.rds.xz

Example
.../by_sample/HC001_B/VDJ_build_HC001_B.rds.xz

3) Patient level VDJ tables
by_patient/<patient_id>/VDJ_build_patient_<patient_id>.rds.xz

Example
.../by_patient/P002/VDJ_build_patient_P002.rds.xz

4) QC summary table
QC_counts_by_sample.tsv

This reports, per sample
capture_id
sample_id
patient_id
tissue
side
n_rows
n_barcodes_demux


Logs

Slurm stdout and stderr
/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/03_logs/


Next step
These patient level VDJ dataframes are ready for downstream steps such as enclone reclonotyping within patient and Af_build in AntibodyForests.
