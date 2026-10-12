# Installation

## Install the pipeline

### Install R dependencies

The TXSeq pipeline requires R>=4.5 and the ``optparse`` ``tximeta`` and ``tximport`` packages.

### Install the Python dependencies

```text
module purge 
unset PYTHONPATH
# Create a python virtual environment Ideally with uv on it 
uv venv txseq_venv --python 3.14
source ./txseq_venv/bin/activate
cd txseq_venv
git clone https://github.com/sansomlab/txseq_snakemake.git
cd txseq_venv
uv pip install -e .
```

## To install the TXSeq R package

Make sure you have the “devtools” and “BiocManager” R libraries pre-installed by running the following command in an R shell:
```text
install.packages(c("devtools","BiocManager"))
```

Install the required R packages by running the following command from the bash shell:
```text
Rscript src/txseq/R/txseq/install.packages.R
```

Install the txseq R library by running the following command from the bash shell:
```text
R CMD INSTALL src/txseq/R/txseq
```
