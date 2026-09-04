"""=========================
Pipeline feature_counts.py
=========================

Overview
--------

This pipeline counts the number of reads mapped to transcript/gene models. It uses the featureCounts algorithm from the `Subread package <https://subread.sourceforge.net>`_ .


Configuration
-------------

The pipeline requires a configured config_feature_counts.yml file.

Inputs
------

The pipeline requires the following inputs

#. samples.tsv: see configuration docs
#. txseq annotations: the location where the ensembl.smk was run to prepare the annotatations.
#. bam files: the location of a folder containing the bam files named by "sample_id".


Requirements
------------

The following software is required:

#. Subread


Output files
------------

The pipeline produces the following outputs:

#. per-sample results: in the "feature.counts.dir" subdirectory
#. A sqlite database: in a file named "csvdb" which contains the per-gene counts.


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
from txseq_snakemake.functions_from_cgat.pipeline.database import concatenate_and_load
from txseq_snakemake.functions_from_cgat.pipeline.parameters import PARAMS
from txseq_snakemake.functions_from_cgat.pipeline.database import to_table
from txseq_snakemake.functions_from_cgat.pipeline.database import load

PARAMS.update(config)

## use txseq tasks to read in sample and library tables
S = T.samples(sample_tsv = config.get("samples"),
        library_tsv = config.get("libraries"))

RESOURCES = config.get("resources", {})

#get bam file name
def bam_file(sample_id):
    bam_fl=os.path.join(config.get("bam_path"),sample_id+".bam")
    return(bam_fl)

#load paired and strand info to make strand param for hisat
paired=config.get("paired")
print(paired)

def featurecounts_strand (sample_id):
    sample = S.samples[sample_id]
    return(sample.featurecounts_strand)


paired=config.get("paired")

if paired:
        paired_options = "-p"
else:
        paired_options = ""

#get options for feature counts
if config["featurecounts"].get("options") is None:
        featurecounts_options = ""
else:
        featurecounts_options = config["featurecounts"].get("options")

# --------------------- Pipeline --------------------------- #
localrules: full
rule full:
    input:
        counts=expand(os.path.join("feature.counts.dir", "{samp}.counts.gz"),samp=S.samples.keys()),
        summary=expand(os.path.join("feature.counts.dir", "{samp}.counts.summary"),samp=S.samples.keys()),
        loadcounts="feature.counts.dir/featurecounts.load",
        gene_counts="feature.counts.dir/featurecounts_counts.tsv.gz",
        load_genecounts="feature.counts.dir/featurecounts_counts.load",
        load_transcript_info=os.path.join("transcript.info.load"),
        n_genes_detected="feature.counts.dir/number.genes.detected.featurecounts.tsv.gz",
        load_n_genes_detected="feature.counts.dir/number.genes.detected.featurecounts.load"

# --------------------- Run feature counts --------------------------- #
rule count:
    input:
        
    output: 
        counts=os.path.join("feature.counts.dir", "{samp}.counts.gz"),
        summary=os.path.join("feature.counts.dir", "{samp}.counts.summary")
    params:
        bam=lambda wildcards: bam_file(wildcards.samp),
        mktemp_template="ctmp.featureCounts.XXXXXXXXXX",
        geneset = os.path.join(config["txseq"].get("annotations"),
                           "api.dir/txseq.geneset.gtf.gz"),
        featurecounts_strand=lambda wildcards: featurecounts_strand(wildcards.samp),
        featurecounts_options=featurecounts_options,
        paired_options=paired_options
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        os.path.join("feature.counts.dir", "{samp}.log")
    shell:
        """
        {{
            counts=`mktemp -p . {params.mktemp_template}`;
                   featureCounts \
                    -a {params.geneset} \
                    -o $counts \
                    -s {params.featurecounts_strand} \
                    -T {threads} \
                    {params.featurecounts_options} \
                    {params.paired_options} \
                    {params.bam};
                    cut -f1,7 $counts \
                    | grep -v "#" | grep -v "Geneid" \
                    | gzip -c > {output.counts};
                    rm $counts;
                    mv ${{counts}}.summary {output.summary};
        }} > {log} 2>&1
        """

# --------------------- Put results in SQL database --------------------------- #
rule loadCounts:
    input:
        expand(os.path.join("feature.counts.dir", "{samp}.counts.gz"),samp=S.samples.keys())
    output:
        "feature.counts.dir/featurecounts.load"
    params:
        job_memory=config["sqlite"].get("himem")
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "feature.counts.dir/featurecounts.load.log"
    run:
            statement = concatenate_and_load(
                    infiles=list(input),
                    outfile=str(output[0]),
                    regex_filename=".*/(.*).counts.gz",
                    has_titles=False,
                    cat="track",
                    header="track,gene_id,counts",
                    job_memory=params.job_memory,
                    options='-i "gene_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

# --------------------- Counts by gene --------------------------- #
rule geneCounts:
    input:
        "feature.counts.dir/featurecounts.load"
    output: 
        "feature.counts.dir/featurecounts_counts.tsv.gz"
    params:
        script=os.path.join(workflow.basedir,
                "python",
                "feature_counts_table.py"),
        database=config["database"].get("db"),
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "feature.counts.dir/featurecounts_counts.log"
    run:
        table = to_table(input[0])

        database = params.database

        statement =f'''
                python {params.script} \
                    --database={params.database} \
                    --table={table} \
                    --outfile={output}
                ''' 
        print("Running statement:\n", statement, flush=True)
        shell(f"{statement} > {log} 2>&1")

# --------------------- Counts by gene into SQL DB --------------------------- #
rule load_gene_counts:
    input:
       "feature.counts.dir/featurecounts_counts.tsv.gz"
    output:
        "feature.counts.dir/featurecounts_counts.load"
    params:
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "feature.counts.dir/featurecounts_counts.load"
    run:
            statement = load(
                    infile=input[0],
                    outfile=str(output[0]),
                    options='-i "gene_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

# --------------------- Transcript info into SQL --------------------------- #
rule load_transcript_info:
    input:
       os.path.join(config["txseq"].get("annotations"),
                          "api.dir/txseq.transcript.info.tsv.gz")
    output:
        os.path.join("transcript.info.load")
    params:
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "transcript.info.load.log"
    run:
            statement = load(
                    infile=input[0],
                    outfile=str(output[0]),
                    job_memory="15G",
                    options='-i "gene_id" -i "transcript_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

# --------------------- Calculate n genes detected --------------------------- #
rule n_genes_detected:
    input:
        "feature.counts.dir/featurecounts.load"
    output: 
        "feature.counts.dir/number.genes.detected.featurecounts.tsv.gz"
    params:
        script=os.path.join(workflow.basedir,
                "python",
                "feature_counts_no_genes_detected.py"),
        database=config["database"].get("db"),
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "feature.counts.dir/number.genes.detected.featurecounts.log"
    run:
        table = to_table(input[0])

        database = params.database

        statement =f'''
                python {params.script} \
                    --database={params.database} \
                    --table={table} \
                    --outfile={output}
                ''' 
        print("Running statement:\n", statement, flush=True)
        shell(f"{statement} > {log} 2>&1")

# --------------------- n genes detected into SQL --------------------------- #
rule loadNGenesDetected:
    input:
       "feature.counts.dir/number.genes.detected.featurecounts.tsv.gz"
    output:
        "feature.counts.dir/number.genes.detected.featurecounts.load"
    params:
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "feature.counts.dir/number.genes.detected.featurecounts.load.log"
    run:
            statement = load(
                    infile=input[0],
                    outfile=str(output[0]),
                    options='-i "sample_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")



