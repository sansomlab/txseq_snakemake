# `Txseq`  ( Snakemake Implementation) 
> **⚠️ WARNING: `txseq`  is currently in pre-release and under active development. 
>Please use it with caution and report any issues by opening a GitHub issue.**
 
Txseq provides a suite of modular workflows for the QC and analysis of bulk RNA-sequencing data on high-performance 
compute clusters. These include pipelines for sanitising Ensembl annotations for RNA-seq analysis, building transcriptome 
indexes and performing quality control with Fastqc and Picard tools. Pipelines for traditional alignment-based mapping with 
HISAT2 and quantitation with `featureCounts` are provided along with a more-modern pseudo-alignment workflow that incorporates 
Salmon and tximeta. Template reports for sequence quality, post-mapping qc analysis, exploratory data analysis and differential 
gene expression are provided as Rmarkdown notebooks.

## Installation

### Install Python Dependencies 

```bash
# Create a python virtual environment Ideally with uv on it 
python -m venv txseq_smk_venv
source ./txseq_smk_venv/bin/activate
git clone https://github.com/sansomlab/txseq_snakemake.git
cd txseq_snakemake
uv pip install -e . 
```
## Using the `txseq` CLI

Installing the package (see above) puts a `txseq` command on your `PATH`. It is
a thin wrapper around Snakemake that locates the packaged workflow files and
launches them on Slurm, so you do not need to reference the `.smk` files by
path.

### Command structure

```bash
txseq <workflow> <task> [options]
```

`<workflow>` is one of the available workflows (listed by `txseq --help`):
`ensembl`, `hisat_index`, `salmon_index`, `fastqc`, `preprocess_quantseq`,
`hisat`, `salmon`, `feature_counts`, and `bamqc`.

`<task>` is either `config` (write out a template configuration file), `full`
(run the complete workflow), or the name of an individual rule to run a single
step.


### Running a workflow

Each workflow is driven by a configuration file that must be present in the
directory from which the command is run, named `<workflow>.yaml`. The typical
sequence is:

```bash
# 1. Move to the directory in which the workflow should run
cd /path/to/run_directory

# 2. Write out the template configuration file (creates ./hisat.yaml)
txseq hisat config

# 3. Edit ./hisat.yaml to set the input paths and parameters for the run

# 4. Run the workflow
txseq hisat full --cores 20 --jobs 25
```

The configuration file is validated before Snakemake is invoked. In
particular, a valid Slurm account must be set under `resources.account`;
the run will not start while this is unset or left as the template
placeholder.


### Options

| Option | Description | Default |
| --- | --- | --- |
| `--cores` | Number of cores available to Snakemake. | `all` |
| `--jobs` | Number of jobs submitted in parallel. | `1` |
| `--dry` | Perform a dry run: construct and display the job graph without executing it. | off |
| `--lock` | Lock the output files on completion by changing their permissions. | off |
| `-v`, `--version` | Print the installed version and exit. | — |

A dry run is recommended before any first execution, to confirm that the
configuration resolves and the expected jobs are planned:

```bash
txseq hisat full --cores 20 --jobs 25 --dry
```
