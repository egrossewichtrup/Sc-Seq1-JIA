# AntibodyForests build and plot pipeline

This repository contains one Slurm ready bash script that builds AntibodyForests objects from per sample Platypus `VDJ_build` outputs, then generates a PDF of lineage tree plots for the largest clonotypes in each sample.

## Script

* `AF_build_and_plot_v1.sh`

## What the script does

### Part 1  Build AntibodyForests objects per sample

For each sample VDJ object found under `INPUT_BY_SAMPLE`, the script:

* Loads `VDJ_build_<sample_id>.rds.xz`
* Selects a clonotype identifier column in this order
  * `clonotype_id_10x_raw`
  * `clonotype_id_10x`
  * `clonotype_id`
  * `clonotype`
* Copies the chosen clonotype column into a working column called `clonotype_id`
* Replaces underscores in clonotype ids with dots to reduce downstream parsing issues
* Chooses a sequence mode based on which trimmed sequence and germline columns exist and contain usable values
  * heavy plus light amino acid trimmed
  * heavy plus light nucleotide trimmed
  * heavy only amino acid trimmed
  * heavy only nucleotide trimmed
* Filters rows to keep only cells with both sequence and germline present for the chosen mode
* Filters clonotypes by number of unique nodes
  * default behavior removes clonotypes that collapse to a single unique node
  * optional behavior keeps those singleton clonotypes
* Calls `AntibodyForests::Af_build()` with `construction.method = "phylo.network.default"`
* Saves one AntibodyForests object per sample and writes a build QC summary

### Part 2  Plot lineage trees per sample

For each saved AntibodyForests object, the script:

* Searches recursively inside the object to collect embedded `igraph` trees
* Parses tree keys to recover sample and clonotype identifiers when possible
* Ranks trees by node count, then by edge count
* Keeps trees with at least `MIN_NODES` nodes
* Plots the top `N_CLONES` trees into a single PDF per sample
* Tries `AntibodyForests::Af_plot_tree()` first
* Falls back to direct `igraph` plotting if `Af_plot_tree()` fails for a given tree
* Writes a tree ranking table per sample plus a plot QC summary

## Inputs

### Required input directory

`INPUT_BY_SAMPLE` must point to a directory containing one subfolder per sample, each holding:

* `VDJ_build_<sample_id>.rds.xz`

By default, the script is configured to use the outputs from your previous step:

* `/hpc/dla_lti/yermanos_group/Sc-Seq1-JIA/2026-02_AntibodyForest_Pipeline/02_results/vdj_build/by_sample`

### Required columns inside each VDJ object

Minimum required:

* `barcode`
* one clonotype column from the list used in Part 1

In addition, at least one usable trimmed sequence and trimmed germline pair must exist:

* `VDJ_sequence_aa_trimmed` with `VDJ_germline_aa_trimmed`, or
* `VDJ_sequence_nt_trimmed` with `VDJ_germline_nt_trimmed`

If light chain columns exist, the script will prefer modes that include both heavy and light:

* `VJ_sequence_aa_trimmed` with `VJ_germline_aa_trimmed`, or
* `VJ_sequence_nt_trimmed` with `VJ_germline_nt_trimmed`

## Software requirements

The script activates a conda environment that must include:

* R
* Platypus
* AntibodyForests
* igraph

The script is designed to run on a Slurm cluster.

## Outputs

All outputs are written under `project_root`, which is set inside the script.

### AntibodyForests objects

`02_results/af_objects/by_sample/<sample_id>/`

* `AntibodyForests_<sample_id>.rds.xz`
* `AntibodyForests_<sample_id>.sessionInfo.txt`
* `AntibodyForests_<sample_id>.run.log`
* `AntibodyForests_<sample_id>.error.txt` (only written on failure)

### Plots

`02_results/plots/by_sample/<sample_id>/`

* `trees_top<N_CLONES>_minNodes<MIN_NODES>_color<COLOR_BY>_label<LABEL_BY>.pdf`
* `clonotype_ranking_<sample_id>.tsv`
* `Af_plot_tree_<sample_id>.run.log`

### QC summaries

`02_results/qc/`

* `Af_build_summary_by_sample.tsv`
* `Af_plot_tree_summary_by_sample.tsv`

### Manifest

`00_config/`

* `samples_for_af_build.tsv`

## How to run

### Run all samples

    sbatch AF_build_and_plot_v1.sh

### Run one sample

Set `ONLY_SAMPLE` to a sample id:

    ONLY_SAMPLE="P001_LK" sbatch AF_build_and_plot_v1.sh

### Overwrite existing AntibodyForests outputs

By default, existing `AntibodyForests_<sample_id>.rds.xz` outputs are skipped. To rebuild:

    FORCE_REBUILD=1 sbatch AF_build_and_plot_v1.sh

### Keep singleton clonotypes

By default, clonotypes that collapse to one unique node are removed. To keep them:

    KEEP_SINGLETONS=1 sbatch AF_build_and_plot_v1.sh

### Change plot parameters

Example:

    N_CLONES=30 MIN_NODES=6 COLOR_BY="tissue" LABEL_BY="name" sbatch AF_build_and_plot_v1.sh

## Configuration variables

You can set these as environment variables at submit time.

Selectors:

* `ONLY_SAMPLE`  empty means run all samples

Build controls:

* `FORCE_REBUILD`  default `0`
* `KEEP_SINGLETONS`  default `0`

Plot controls:

* `N_CLONES`  default `20`
* `MIN_NODES`  default `4`
* `COLOR_BY`  default `isotype`
* `LABEL_BY`  default `size`
* `EDGE_LABEL`  default `none`  set to `original` to attempt edge labels from graph attributes
* `NODE_SIZE`  default `expansion`  set to `fixed` for constant sizing in the fallback plot
* `SHOW_INNER_NODES`  default `FALSE`

Note: Supported values for `COLOR_BY` and related options depend on what node metadata exist in your VDJ objects and on the function signature of your installed AntibodyForests version.

## Robustness details

* The script checks the formal arguments of `Af_plot_tree()` and only passes arguments that exist in your installed AntibodyForests version.
* If `Af_plot_tree()` cannot plot a specific tree, the script will still plot it using a direct `igraph` fallback.
* Tree ranking is based on node count first, then edge count, after filtering by `MIN_NODES`.

## Troubleshooting

If the build step fails for a sample:

* Check `02_results/af_objects/by_sample/<sample_id>/AntibodyForests_<sample_id>.run.log`
* If present, check `AntibodyForests_<sample_id>.error.txt`
* Common causes are missing `barcode`, missing clonotype columns, or missing trimmed sequence and germline columns

If the plot step produces no trees:

* Check `02_results/plots/by_sample/<sample_id>/clonotype_ranking_<sample_id>.tsv`
* Check `02_results/qc/Af_plot_tree_summary_by_sample.tsv`
* A common reason is that no trees pass `MIN_NODES`

If the PDF exists but looks cluttered:

* Try `LABEL_BY="none"`
* Try `COLOR_BY` set to a column that exists for nodes, such as `isotype` or `tissue`

## Reproducibility

Each sample folder includes a saved `sessionInfo()` dump from the build step so you can track the exact R package versions used to construct the AntibodyForests objects.