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


#rule all
rule full:
    input:
        expand(os.path.join(config["rundir"], "hisat.dir","first.pass.dir", "{samp}.novel.splice.sites.txt.gz"), samp=S.samples.keys()),


## first pass of hisat. 
rule first_pass:
    input:
        os.path.join(config["fq_dir"],"{samp}_1.fq.gz")
    output:
        os.path.join(config["rundir"], "hisat.dir","first.pass.dir", "{samp}.novel.splice.sites.txt.gz")
    params:
        index=config["txseq"].get("hisat_index"),
    #    threads=config["hisat_thread"],
        opts=config["hisat"].get("options"),
        strand_param=strand_param,
        fastq_input=lambda wildcards: fastq_input(wildcards.samp),
        known_ss=ss_file,
    shell:
        """
        hisat2 \
                    -x {params.index} \
                    {params.fastq_input} \
                    --known-splicesite-infile {params.known_ss} \
                    --novel-splicesite-outfile {output} \
                    {params.strand_param} \
                    {params.opts} \
                    -S /dev/null; 
                gzip {output}
        """