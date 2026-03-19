# Cambodia BaseSpace example

This is intended as the repo's quick start, so we provide a little extra documentation here in case you got stuck on the main quick start. As always, make sure you have cloned the repository, installed the environment, and activated the environment with:

```
conda activate mlip
```

A file of IDs has been provided for you, though it's a good exercise to look at the data directory and try to create this yourself.

An example metadata spreadsheet that will run has also been provided for you. Again, try to get to this point yourself, but check if you feel stuck.

Run preprocess to initialize your analysis:

```
python mlip/dataflow.py preprocess -f examples/cambodia-basespace/ids.txt --analysis cambodia
```

Edit `cambodia/config.yml` to set `reference: "h5n1"` and `data_root_directory: "examples/cambodia-basespace"`. Populate sample IDs and replicates in `cambodia/metadata.tsv`. Then run:

```
python mlip/dataflow.py configure
```

Run the pipeline:

```
snakemake -j $NUMBER_OF_JOBS all
```

`$NUMBER_OF_JOBS` should be at least 1, and no more than the number of cores on your computer.

Run to explore outputs:

```
python mlip/visualization.py
```
