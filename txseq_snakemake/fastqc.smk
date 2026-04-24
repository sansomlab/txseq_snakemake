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

## use txseq tasks to read in sample and library tables
S = T.samples(sample_tsv = config.get("samples"),
        library_tsv = config.get("libraries"))


def fastq_path(samp): 
    return(S.fastqs[samp]['fastq_path'])

def make_fastq_outpath(samp):
    path=os.path.join("fastqc.dir", samp)
    if not os.path.exists(path):
        os.makedirs(path)
    return(path)


contaminants = ""
if config["fastqc"].get("contaminants") != "default":
        contaminants = "--contaminants " + config["fastqc"].get("contaminants")
    
adaptors = ""
if config["fastqc"].get("adaptors") != "default":
        adaptors = "-a " + config["fastqc"].get("adaptors")
        
limits = ""
if config["fastqc"].get("limits") != "default":
    limits = "-l " + config["fastqc"].get("limits")

rule full:
    input:
        expand(os.path.join("fastqc.dir", "{samp}.sentinel"), samp=S.fastqs.keys())

rule fastqc:
    input:
        
    output: 
        os.path.join("fastqc.dir", "{samp}.sentinel")
    params:
        fastq_path=lambda wildcards: fastq_path(wildcards.samp),
        out_path = lambda wildcards: make_fastq_outpath(wildcards.samp),
        adaptors=adaptors,
        limits=limits,
        contaminants=contaminants,
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        os.path.join("fastqc.dir", "{samp}.log")
    shell:
        """
        {{
            fastqc \
                   -o {params.out_path} \
                   --extract \
                   {params.contaminants} \
                   {params.adaptors} \
                   {params.limits} \
                   {params.fastq_path};
            touch {output}
        }} > {log} 2>&1
        """

