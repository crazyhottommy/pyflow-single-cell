# snakemake pipeline for processing multi-sample scATACseq data


```bash

git clone https://github.com/crazyhottommy/pyflow-single-cell
cd scATAC

# activate snakemake env
conda activate MAESTRO

## generate a dummy samples.json file
../utils/sample2json.py scatac ./data/

```

A `samples.json` file will be generated.  


The `sample2json.py` script uses regular expression to extract file names. This requires
the files are named consistently. This is fine for 10x data downloaded from their website. In the future, a `meta.txt` file need to be used if the users do not have consistent file naming.

```bash

## dry run
snakemake -np

```