# Bioinformatics pipeline for Illumina viral deep sequencing

WARNING: This repository is a work in progress.

This `README.md` is intended to be a quickstart overview. For a deeper understanding of this pipeline, please see [our full documentation](./DOCUMENTATION.md).

## Installation

Requires [Bioconda](https://bioconda.github.io/) and [Git](https://git-scm.com/). We recommend [Miniconda](https://docs.anaconda.com/miniconda/) be used as your conda distribution.

Create an environment with the tools used by this pipeline:
```
conda create -n mlip python=3.12 pandas=2 altair biopython bedtools bcftools bowtie2 multiqc samtools trimmomatic snakemake=8.27 varscan entrez-direct seqkit sed csvkit perbase vapor mafft ivar
```

## Usage

Usage instructions assume that you've successfully followed the [installation instructions](#installation), and read about and adhere to the [conventions](#conventions) used by this software. Further, it assumes a basic understanding of command line interfaces, as well as [conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/en/stable/).

If it's your first run, we've already prepared data and encourage you to use our [Cambodia BaseSpace example](./examples/cambodia-basespace). We have also prepared other [examples](./examples).

### Quick start
Suppose you have several FASTQs downloaded to a folder that you'd like to analyze for a project called `MyAnalysis`.

Get a copy of the code and set up your environment for analysis:

```
git clone https://github.com/moncla-lab/illumina-pipeline MyAnalysis
cd MyAnalysis
conda activate mlip
```

There are three main concerns when configuring the pipeline to run:

- copying and editing the **configuration** template, and pointing it towards existing data 
- choosing a **reference**, either in the configuration, in the references table, or by using a custom one
- generating and filling in an appropriate **metadata** spreadsheet

The user is encouraged to **repeatedly** run

```
python mlip/dataflow.py check
```

to receive feedback on where they are in the configuration process.

#### Configuration file

Copy the file `config.yml.template` to `config.yml`. Make appropriate edits, which will likely involve adjusting only `reference` (what's used as the reference sequence) and the `data_root_directory` (where data was downloaded) for a first run.

#### References
We have [predefined references](./references.tsv). The simplest use case is to choose a key from the reference column to populate the config. The user can override these by defining their own with segments pulled from Genbank or using a custom reference. There is [extended documentation on references](./DOCUMENTATION.md#references) for more detail.

#### Metadata

To move data out of the data folder and into this pipeline, a text file of sequencing experiment IDs corresponding to each FASTQ dataset, one per line, must be created at a known path we'll call `/path/to/fastqDatasetIDs.txt`. Note we say FASTQ dataset as this is technically a pair of FASTQ files. It's best to pull these from either BaseSpace sample sheet CSVs or SRA accessions.

The first **preprocess** step looks at these IDs and builds a metadata spreadsheet at `data/metadata.tsv`:
```
python mlip/dataflow.py preprocess -f /path/to/sequencingExperimentIDs.txt
```

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

For emphasis, the `check` step above will attempt to give feedback at any step of the configuration.

Once the metadata spreadsheet is fully populated, the following command moves data into the pipeline in an organized manner for further analysis.

```
python mlip/dataflow.py flow
```

Note that the default mode for the above command is to handle BaseSpace data, and there is an SRA mode that uses the `--sra-mode` flag.

### Run the pipeline

With data situated, the pipeline can be ran as:
```
snakemake -j $NUMBER_OF_JOBS all
```

$NUMBER_OF_JOBS should be at least 1, and no more than the number of cores on your computer. After this, the `data` directory should be filled with lots of files of various formats, many which contain relevant virological information.

### Outputs

To bring up a directory tree of the `data` directory where you will find files of interest and be able to view certain plots, run:

```
python mlip/visualization.py
```

Alternatively, just explore the data directory from your desktop. All paths below are assumed to be relative to the `data` directory. Anything enclosed in brackets are [Snakemake wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-wildcards) with further explanation in documentation. Relevant outputs include:

| File description                        | File path                                                     |
| --------------------------------------- | ------------------------------------------------------------- |
| Consensus sequences for a given segment | `{segment}.fasta`                                             |
| Protein sequences for a given gene      | `protein/{gene}.fasta`                                        |
| Annotated, merged variants              | `variants.tsv`                                                |
| Project wide overview of coverage       | `coverage-report.tsv`                                         |
| Zip of all small files                  | `project.zip`                                                 |
| Plot of intrahost variants for a sample | `{sample}/ml.html`                                            |
| Replicate, mapping specific coverage    | `{sample}/replicate-{replicate}/{mapping_stage}/coverage.tsv` |

We also have a [more comprehensive list of outputs](DOCUMENTATION.md#5-output-produced).

## Conventions

This pipeline makes some assumptions about how data is organized. In particular, BaseSpace/sequencing experiment IDs are of the following format:

```
{SAMPLE}_Seq{#}
```

See the example metadata above. This helps automatically populate the sample associated to a sequencing experiment, and helps keep the data in BaseSpace and in the pipeline in sync.

There are analogous formats for SRA, and we expect users to put their data in one format or another to work with this pipeline.

This gives an overall gist of the pipeline. For further explanation, please consult [our documentation](./DOCUMENTATION.md).
