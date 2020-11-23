# extreme-memorizers

Collection of snakemake pipeline and jupyter notebooks used for the analysis
of the extreme-memorizers carriers against a set of WGS data.

## Pipeline

This repository contains a set of snakemake scripts that can be used to
generate the output data. The `results.ipynb` playbook contains some examples
for working with the generated variant tables.

Make sure you have snakemake and (mini)conda installed to run this pipeline.

## Executing

To execute the pipeline, please adjust the configuratoin file paths to their
corresponding values. The pipeline will search all `input_dirs` for files
matching the given glob patterns (`search-patterns`) and annotate and filter
all these scripts.