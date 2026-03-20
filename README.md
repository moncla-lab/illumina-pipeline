# Bioinformatics pipeline for Illumina viral deep sequencing

WARNING: This repository is a work in progress.

This `README.md` is intended to be a quickstart overview. For a deeper understanding of this pipeline, please see [our full documentation](./DOCUMENTATION.md).

## Installation

Requires [Bioconda](https://bioconda.github.io/) and [Git](https://git-scm.com/). We recommend [Miniconda](https://docs.anaconda.com/miniconda/) be used as your conda distribution.

Create an environment with the tools used by this pipeline:

```
conda env create -f environment.yml
```

## Usage

Usage instructions assume that you've successfully followed the [installation instructions](#installation), and read about and adhere to the [conventions](#conventions) used by this software. Further, it assumes a basic understanding of command line interfaces, as well as [conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/).

If it's your first run, we've already prepared data and encourage you to use our [Cambodia BaseSpace example](./examples/cambodia-basespace). We have also prepared other [examples](./examples).

### Quick start
Suppose you have several FASTQs downloaded to a folder that you'd like to analyze for a project called `MyAnalysis`.

Get a copy of the code and set up your environment for analysis:

```
git clone https://github.com/moncla-lab/illumina-pipeline
cd illumina-pipeline
conda activate mlip
```

There are three main concerns when configuring the pipeline to run:

- initializing an analysis with a **configuration** file
- choosing a **reference**, either in the configuration, in the references table, or by using a custom one
- generating and filling in an appropriate **metadata** spreadsheet

#### Initialize your analysis

First, create a text file with your sequencing experiment IDs (one per line). Then run **preprocess** to initialize your analysis:

```
python mlip/dataflow.py preprocess -f /path/to/sequencingExperimentIDs.txt -a MyAnalysis
```

This creates:

- `MyAnalysis/config.yml` - Configuration file for this analysis
- `MyAnalysis/metadata.tsv` - Metadata spreadsheet to edit

#### Configuration file

Edit `MyAnalysis/config.yml` to set the required values:

- `reference`: Reference key from `references.tsv` or path to custom reference ZIP (required)
- `data_root_directory`: Path where you downloaded data from the BaseSpace downloader (required)
- Coverage thresholds if needed

For a complete list of configuration options, see [Configuration documentation](./DOCUMENTATION.md#2-configuration).

#### References
We have [predefined references](./references.tsv). The simplest use case is to choose a key from the reference column to populate the config. The user can override these by defining their own with segments pulled from Genbank or using a custom reference. There is [extended documentation on references](./DOCUMENTATION.md#references) for more detail.

#### Metadata

The **preprocess** step above created a metadata spreadsheet at `MyAnalysis/metadata.tsv`.

The pipeline will do its best to assign sample IDs to each sequencing experiment ID. The user should open the metadata file above, inspect that these sample IDs were correctly assigned, assign any that are missing, and assign a replicate to each sequencing experiment ID.

For clarity, the input at this step may look something like:

```
bv_w1_Seq1
bv_w1_Seq2
rf_Seq1
rf_Seq2
rf_Seq3
```

while a completed metadata sheet will look like:

| SequencingId | SampleId | Replicate |
| ------------ | -------- | --------- |
| bv\_w1_Seq1  | bv_w1    | 1         |
| bv\_w1_Seq2  | bv_w1    | 2         |
| rf_Seq1      | rf       | 1         |
| rf_Seq2      | rf       | 2         |
| rf_Seq3      | rf       | 1         |

#### Configure and validate

Once the metadata spreadsheet is fully populated and you've edited the config, run **configure** to validate your setup and prepare data for the pipeline:

```
python mlip/dataflow.py configure
```

This command:

1. Validates your configuration and metadata
2. If all checks pass, automatically prepares data for the pipeline
3. Reports any issues that need to be fixed

By default, `configure` uses the most recently modified `config.yml`. To target a specific analysis, use `-a`:

```
python mlip/dataflow.py configure -a MyAnalysis
```

Note: For SRA data instead of BaseSpace, use `-s` / `--sra-mode` flag.

### Run the pipeline

With data situated, the pipeline can be ran as:
```
snakemake -j $NUMBER_OF_JOBS all
```

`$NUMBER\_OF\_JOBS` should be at least 1, and no more than the number of cores on your computer. After this, your analysis output directory (named per the `analysis` parameter in `config.yml`) should be filled with lots of files of various formats, many which contain relevant virological information.

**Important**: Always review the consensus remapping differences report (`consensus_summary_report.tsv`) after running the pipeline to identify positions where consensus sequences changed between remapping iterations, which indicates potential mapping inconsistencies that should be investigated. When this occurs, consult the [additional documentation](DOCUMENTATION.md#check_consensus_summary).

### Outputs

To bring up a directory tree of your analysis output directory where you will find files of interest and be able to view certain plots, run:

```
python mlip/visualization.py
```

Alternatively, just explore the analysis directory from your desktop. All paths below are assumed to be relative to your analysis output directory (e.g., `stephen-test/` if `analysis: "stephen-test"` in `config.yml`). Anything enclosed in brackets are [Snakemake wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-wildcards) with further explanation in documentation. Relevant outputs include:

| File description                        | File path                                                     |
| --------------------------------------- | ------------------------------------------------------------- |
| Consensus sequences for a given segment | `{segment}.fasta`                                             |
| Protein sequences for a given gene      | `protein/{gene}.fasta`                                        |
| Annotated, merged variants              | `variants.tsv`                                                |
| Project wide overview of coverage       | `coverage-report.tsv`                                         |
| Consensus remapping differences report  | `consensus_summary_report.tsv`                                |
| Zip of all small files                  | `project.zip`                                                 |
| Plot of intrahost variants for a sample | `{sample}/ml.html`                                            |
| Replicate, mapping specific coverage    | `{sample}/replicate-{replicate}/{mapping_stage}/coverage.tsv` |
| Final BAM files (preserved from last remapping) | `{sample}/replicate-{replicate}/final.bam`               |

We also have a [more comprehensive list of outputs](DOCUMENTATION.md#5-output-produced).

## Conventions

This pipeline makes some assumptions about how data is organized. In particular, BaseSpace/sequencing experiment IDs are of the following format:

```
{SAMPLE}_Seq{#}
```

See the example metadata above. This helps automatically populate the sample associated to a sequencing experiment, and helps keep the data in BaseSpace and in the pipeline in sync.

There are analogous formats for SRA, and we expect users to put their data in one format or another to work with this pipeline.

This gives an overall gist of the pipeline. For further explanation, please consult [our documentation](./DOCUMENTATION.md).
