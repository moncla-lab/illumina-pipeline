# Cambodia BaseSpace example

This is intended as the repo's quick start, so we provide a little extra documentation here in case you got stuck on the main quick start. After installing, edit the config to contain:

```
reference: "h5n1"
data_root_directory: "~/path/to/illumina-pipeline/examples/cambodia-basespace"
```

A file of IDs has been provided for you, though it's a good exercise to look at the data directory and try to this create yourself.

An example metadata spreadsheet that will run has also been provided for you. Again, try to get to this point yourself, but check if you feel stuck.

Run

```
python mlip/dataflow preprocess -f ~/path/to/cambodia-basespace/ids.txt
```

to generate the metadata. Populated sample IDs and replicates in the metadata. Run

```
python mlip/dataflow flow
```

to show the repository where the data lies. Run

```
snakemake -j $NUMBER_OF_JOBS all
```

to extract data.

Run

```
python mlip/visualization.py
```

to explore outputs.