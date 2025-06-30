# Cambodia BaseSpace example

This is intended as the repo's quick start, so we provide a little extra documentation here in case you got stuck on the main quick start. As always, make sure you have cloned the repository, installed the environment, and activated the environment with:

```
conda activate mlip
```

After installing, edit the config to contain:

```
reference: "h5n1"
data_root_directory: "~/path/to/illumina-pipeline/examples/cambodia-basespace"
```

A file of IDs has been provided for you, though it's a good exercise to look at the data directory and try to this create yourself.

An example metadata spreadsheet that will run has also been provided for you. Again, try to get to this point yourself, but check if you feel stuck.

Run

```
python mlip/dataflow.py preprocess -f ~/path/to/cambodia-basespace/ids.txt
```

to generate the metadata. This will create a metadata spreadsheet at `./data/metadata.tsv`. Populate sample IDs and replicates in the metadata. Run

```
python mlip/dataflow.py flow
```

to show the repository where the data lies. Run

```
snakemake -j $NUMBER_OF_JOBS all
```

to extract data. `$NUMBER_OF_JOBS` should be at least 1, and no more than the number of cores on your computer.

Run

```
python mlip/visualization.py
```

to explore outputs.
