"""=======================
pipeline_hisat_index.py
=======================

Overview
--------

This pipeline makes a hisat2 HGFM index with transcripts.

Configuration
-------------

The pipeline requires a configured config_hisat_index.yml file.



Input files
-----------

For building the  index the pipeline requires:

#. A gzip compressed FASTA file containing the genome primary assembly sequences.
#. A gzip compressed GTF file containing the transcript information.

The location of these files must be specified in the 'config_hisat_index.yml' file.

Requirements
------------

The following software is required:

#. Hisat2


Output files
------------

The pipeline generates a hisat2 HGFM index in the current folder.

Code
====

"""

import sys
import shutil
import os
from pathlib import Path
import glob
import sqlite3

import pandas as pd
import numpy as np

# import local pipeline utility functions
import txseq_snakemake.tasks as T

RESOURCES = config.get("resources", {})

# --------------------- Pipeline --------------------------- #
localrules: full
#rule all
rule full:
    input:
        os.path.join(config["rundir"],"genome_tran.sentinel")

# --------------------- Identify splice sites --------------------------- #
## splice sites. 
rule splice_sites:
    input:
        config["geneset"]
    output:
        os.path.join(config["rundir"],"genome.ss")
    params:
    threads: 1
    resources:
        **RESOURCES,
    log:
        "splice_sites.log"
    shell:
        """
        {{
            hisat2_extract_splice_sites.py \
            <(zcat {input}) \
            > {output} 
        }} > {log} 2>&1
        """

# --------------------- Identify exons --------------------------- #
## exons 
rule exons:
    input:
        config["geneset"]
    output:
        os.path.join(config["rundir"],"genome.exon")
    params:
    threads: 1
    resources:
        **RESOURCES,
    log:
        "exons.log"
    shell:
        """
        {{
            hisat2_extract_exons.py  \
            <(zcat {input}) \
            > {output} 
        }} > {log} 2>&1
        """

# --------------------- Build trnascriptome index --------------------------- #
## transcriptome index
rule transcriptome_index:
    input:
        ss=os.path.join(config["rundir"],"genome.ss"),
        exon=os.path.join(config["rundir"],"genome.exon")
    output:
        os.path.join(config["rundir"],"genome_tran.sentinel")
    params:
        primary_assembly=config["primary_assembly"],
        outs=os.path.join(config["rundir"],"genome_tran"), 
    threads: 
        config["index_threads"]
    resources:
        mem_mb=config["index_memory"],
        runtime=config["index_runtime"],
    log:
        "transcriptome_index.log"
    shell:
        """
        {{
            zcat {params.primary_assembly} > genome.fa;
                   hisat2-build \
                   -p {threads} \
                   --exon {input.exon} \
                   --ss {input.ss} \
                   genome.fa \
                   genome_tran; 
                   rm genome.fa;
            touch {output}
        }} > {log} 2>&1
        """

