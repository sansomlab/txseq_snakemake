"""=================
Pipeline hisat.py
=================

Overview
--------

This pipeline quantifies gene expression from FASTQ files using `Hisat2 <http://daehwankimlab.github.io/hisat2/>`_. 


Configuration
-------------

The pipeline requires a configured :file:`pipeline_hisat.yml` file.


Inputs
------

The pipeline requires the following inputs

#. samples.tsv: see :doc:`Configuration files<configuration>`
#. libraries.tsv: see :doc: `Configuration files<configuration>`
#. txseq annotations: the location where the :doc:`pipeline_ensembl.py </pipelines/pipeline_ensembl>` was run to prepare the annotatations.
#. Hisat index: a hisat2 index built with :doc:`pipeline_hisat_index.py </pipelines/pipeline_hisat_index>`.


Requirements
------------

The following software is required:

#. Hisat2


Output files
------------

The pipeline produces the following outputs:

#. bam files: these are found in the hisat.dir sub-folder and are named by sample_id.


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

##function- takes sample id and puts the strand input together as needed for hisat
def fastq_input(sample_id):
    sample = S.samples[sample_id]
    if sample.paired:
        fastq_input = "-1 " + ",".join(sample.fastq["read1"]) +\
                      " -2 " + ",".join(sample.fastq["read2"])
    else:
        fastq_input = "-U " + ",".join(sample.fastq["read1"])
    return(fastq_input)

## use txseq tasks to read in sample and library tables
S = T.samples(sample_tsv = config.get("samples"),
                 library_tsv = config.get("libraries"))
print(S)
#load paired and strand info to make strand param for hisat
paired=config.get("paired")
print(paired)
strand=config.get("strand")
print(strand)


if strand=="forward":
    if paired:
        strand_param = "--rna-strandness FR"
    else:
        strand_param = "--rna-strandness F"
elif strand=="reverse":
    if paired:
        strand_param = "--rna-strandness RF"
    else:
        strand_param = "--rna-strandness R"
else:
    strand_param = ""
print(strand_param)


#known splice sites file
ss_file = os.path.join(os.path.dirname(config["txseq"].get("hisat_index")),
                               "genome.ss")
print("ss file is "+ss_file)

RESOURCES = config.get("resources", {})


# --------------------- Pipeline --------------------------- #
#rule all
rule full:
    input:
        expand(os.path.join("hisat.dir", "{samp}.bam"), samp=S.samples.keys())

# --------------------- First pass --------------------------- #
## first pass of hisat. 
rule first_pass:
    input:
        os.path.join(config["fq_dir"],"{samp}_1.fq.gz")
    output:
        os.path.join("hisat.dir","first.pass.dir", "{samp}.novel.splice.sites.txt.gz")
    params:
        index=config["txseq"].get("hisat_index"),
    #    threads=config["hisat_thread"],
        opts=config["hisat"].get("options"),
        strand_param=strand_param,
        fastq_input=lambda wildcards: fastq_input(wildcards.samp),
        known_ss=ss_file,
        novel_ss_outfile=os.path.join("hisat.dir","first.pass.dir", "{samp}.novel.splice.sites.txt"),
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        os.path.join("hisat.dir","first.pass.dir", "{samp}.log")
    shell:
        """
        {{
        hisat2 \
                    -x {params.index} \
                    --threads {threads} \
                    {params.fastq_input} \
                    --known-splicesite-infile {params.known_ss} \
                    --novel-splicesite-outfile {params.novel_ss_outfile} \
                    {params.strand_param} \
                    {params.opts} \
                    -S /dev/null; 
                gzip {params.novel_ss_outfile} 
        }} > {log} 2>&1
        """

# --------------------- Collect splice sites --------------------------- #
## collect all novel splice sites files into one
rule collect_sites:
    input:
        expand(os.path.join("hisat.dir","first.pass.dir", "{samp}.novel.splice.sites.txt.gz"), samp=S.samples.keys()),
    output:
        os.path.join("hisat.dir", "annotations","novel.splice.sites.txt")
    threads: 1
    resources:
        cpus_per_task=1
    log:
        os.path.join("hisat.dir", "annotations","novel.splice.sites.log")
    shell:
        """
        {{
        mkdir -p tmp.dir;
        sort_dir=`mktemp -d -p tmp.dir`;
        zcat {input} \
        | sort -k1,1 -T $sort_dir \
        | uniq \
        > {output};
        rm -rf $sort_dir
        }} > {log} 2>&1
        """

# --------------------- Second pass --------------------------- #
## second pass of hisat. 
rule second_pass:
    input:
        sample_int=os.path.join(config["fq_dir"],"{samp}_1.fq.gz"),
        splice_file=os.path.join("hisat.dir", "annotations","novel.splice.sites.txt")
    output:
        os.path.join("hisat.dir", "{samp}.bam")
    params:
        index=config["txseq"].get("hisat_index"),
    #    threads=config["hisat_thread"],
        opts=config["hisat"].get("options"),
        strand_param=strand_param,
        fastq_input=lambda wildcards: fastq_input(wildcards.samp),
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        os.path.join("hisat.dir", "{samp}.log")
    shell:
        """
        {{
        mkdir -p tmp.dir;
            sort_dir=`mktemp -d -p tmp.dir`;
                hisat2 \
                        -x {params.index} \
                        --threads {threads} \
                        {params.fastq_input} \
                        --novel-splicesite-infile {input.splice_file} \
                        {params.strand_param} \
                        {params.opts} \
                    |    samtools view - -bS \
                    | samtools sort - -T $sort_dir -o {output};
                    samtools index {output};
                    rm -rf $sort_dir;
        }} > {log} 2>&1
        """
