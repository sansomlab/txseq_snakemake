"""==================
Pipeline salmon.py
==================

Overview
--------

This pipeline quantifies gene expression from FASTQ files using `Salmon <https://github.com/COMBINE-lab/salmon>`_. 


Configuration
-------------

The pipeline requires a configured :file:`config_salmon.yml` file.


Inputs
------

The pipeline requires the following inputs

#. samples.tsv: see documentation configuration page
#. libraries.tsv: see see documentation configuration page
#. txseq annotations: the location where the ensembl.smk pipeline was run to prepare the annotatations.
#. Salmon index: the location of a salmon index built with salmon_index.smk.


Requirements
------------

The following software is required:

#. Salmon

Output files
------------

The pipeline produces the following outputs:

#. per-sample salmon quantification results in the "salmon.dir" folder.
#. a csvdb sqlite database that contains tables of gene and transcript counts and TPMs.
#. a gene level summary in a genesummarizedexperiment in the "tximeta" folder. 
#. information on the number of genes detected in each sample.

.. note::

    It is strongly recommended to parse the raw Salmon results using the `tximport <https://bioconductor.org/packages/release/bioc/html/tximport.html>`_ Bioconductor R package for downstream analysis.
    

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

#General pipeline prep

PARAMS.update(config)
RESOURCES = config.get("resources", {})

## use txseq tasks to read in sample and library tables
S = T.samples(sample_tsv = config.get("samples"),
                 library_tsv = config.get("libraries"))

#function to get fq input files
def fastq_input(sample_id):
    sample = S.samples[sample_id]
    if sample.paired:
        fastq_input = "-1 " + ",".join(sample.fastq["read1"]) +\
                      " -2 " + ",".join(sample.fastq["read2"])
    else:
        fastq_input = "-r " + ",".join(sample.fastq["read1"])
    return(fastq_input)


#get salmon quant options
salmon_quant_options = ''
if not config.get('salmon_quant_options') is None:
        salmon_quant_options = config.get('salmon_quant_options')

#get salmon library type
def libtype(sample_id):
    sample = S.samples[sample_id]
    return(sample.salmon_libtype)

# --------------------- Pipeline --------------------------- #
localrules: full
rule full:
    input:
        quant=expand(os.path.join("salmon.dir","{samp}","quant.sf"), samp=S.samples.keys()),
        load_transcripts=os.path.join("salmon.dir", "salmon.transcripts.load"),
        gene_quant=os.path.join("salmon.dir", "salmon.genes.load"),
        salmon_genes_tpms=os.path.join("salmon.dir", "salmon.genes.tpms.txt.gz"),
        salmon_transcripts_tpms=os.path.join("salmon.dir", "salmon.transcripts.tpms.txt.gz"),
        salmon_transcripts_tpms_load=os.path.join("salmon.dir", "salmon.transcripts.tpms.load"),
        salmon_genes_tpms_load=os.path.join("salmon.dir", "salmon.genes.tpms.load"),
        tximeta=os.path.join("tximeta.dir","tximeta.RDS"),
        transcript_infoload=os.path.join("transcript.info.load"),
        ngenes_detected=os.path.join("qc.dir", "number.genes.detected.salmon.tsv.gz"),
        load_ngenes=os.path.join("qc.dir","qc_no_genes_salmon.load")

# --------------------- Run salmon quant --------------------------- #
rule quant:
    input:
    output:
        quant=os.path.join("salmon.dir","{samp}","quant.sf"),
        quant_genes=os.path.join("salmon.dir","{samp}","quant.genes.sf")
    params:
        fastq_input=lambda wildcards: fastq_input(wildcards.samp),
        salmon_quant_options=salmon_quant_options,
        libtype=lambda wildcards: libtype(wildcards.samp),
        tx2gene = os.path.join(config["txseq"].get("annotations"),"api.dir/txseq.transcript.to.gene.map"),
        txseq_salmon_index=config["txseq"].get("salmon_index"),
        output_dir=os.path.join("salmon.dir","{samp}")
    threads:
        config["salmon"].get("threads"),
    resources:
        mem_mb=config["salmon"].get("memory"),
        runtime=config["salmon"].get("runtime"),
    log:
                os.path.join("salmon.dir","{samp}_quant.log")
    shell:
        """
        {{
            salmon quant -i {params.txseq_salmon_index} \
                                -p {threads} \
                                -g {params.tx2gene} \
                                {params.salmon_quant_options} \
                                -l {params.libtype} \
                                {params.fastq_input} \
                                -o {params.output_dir}
        }} > {log} 2>&1
        """

# --------------------- Load into SQL --------------------------- #
rule loadSalmonTranscriptQuant:
    input:
       quant=expand(os.path.join("salmon.dir","{samp}","quant.sf"), samp=S.samples.keys()),
    output:
        os.path.join("salmon.dir", "salmon.transcripts.load")
    params:
        job_mem=config["sqlite"].get("himem")
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("salmon.dir", "salmon.transcripts.load.log")
    run:
            statement = concatenate_and_load(
                    infiles=list(input),
                    outfile=str(output[0]),
                    regex_filename=".*/(.*)/quant.sf",
                    cat="sample_id",
                    job_memory={params.job_mem},
                    options="-i Name -i sample_id"
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")


rule loadSalmonGeneQuant:
    input:
       quant=expand(os.path.join("salmon.dir","{samp}","quant.genes.sf"), samp=S.samples.keys()),
       preload1=os.path.join("salmon.dir", "salmon.transcripts.load")
    output:
        os.path.join("salmon.dir", "salmon.genes.load")
    params:
        job_mem=config["sqlite"].get("himem")
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("salmon.dir", "salmon.genes.load.log")
    run:
            statement = concatenate_and_load(
                    infiles=list(input.quant),
                    outfile=str(output[0]),
                    regex_filename=".*/(.*)/quant.genes.sf",
                    cat="sample_id",
                    job_memory={params.job_mem},
                    options="-i Name -i sample_id"
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

# --------------------- Get TPMs --------------------------- #
rule salmon_genes_tpms:
    input:
       os.path.join("salmon.dir", "salmon.genes.load")
    output:
        os.path.join("salmon.dir", "salmon.genes.tpms.txt.gz")
    params:
        database=config["database"].get("db"),
        script=os.path.join(workflow.basedir,
                "python",
                "salmon_fetch_tpms.py"),
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("salmon.dir", "salmon.genes.tpms.load.log")
    run:
            table = to_table(input[0])
            outfile=os.path.join("salmon.dir", "salmon.genes.tpms.txt.gz")

            statement = f"""
                python {params.script} \
                   --database={params.database} \
                   --table={table} \
                   --idname=gene_id \
                   --outfile={outfile}
                """

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

rule salmon_transcripts_tpms:
    input:
       os.path.join("salmon.dir", "salmon.transcripts.load")
    output:
        os.path.join("salmon.dir", "salmon.transcripts.tpms.txt.gz")
    params:
        database=config["database"].get("db"),
        script=os.path.join(workflow.basedir,
                "python",
                "salmon_fetch_tpms.py"),
    threads: 1
    resources:
        mem_mb=config["salmon"].get("memory"),
        runtime=config["salmon"].get("runtime"),
    log:
        os.path.join("salmon.dir", "salmon.transcripts.tpms.load.log")
    run:
            table = to_table(input[0])
            outfile=os.path.join("salmon.dir", "salmon.transcripts.tpms.txt.gz")

            statement = f"""
                python {params.script} \
                   --database={params.database} \
                   --table={table} \
                   --idname=transcript_id \
                   --outfile={outfile}
                """

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

# --------------------- Load TPMs into SQL --------------------------- #
rule loadSalmonTranscriptTPMs:
    input:
       tpms=os.path.join("salmon.dir", "salmon.transcripts.tpms.txt.gz"),
       preload1=os.path.join("salmon.dir", "salmon.genes.load")
    output:
        os.path.join("salmon.dir", "salmon.transcripts.tpms.load")
    params:
        job_mem=config["sqlite"].get("himem")
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("salmon.dir", "salmon.transcripts.tpms.load.log")
    run:
            statement = load(
                    infile=input.tpms,
                    outfile=str(output[0]),
                    job_memory={params.job_mem},
                    options="-i transcript_id"
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

rule loadSalmonGeneTPMs:
    input:
       tpms=os.path.join("salmon.dir", "salmon.genes.tpms.txt.gz"),
       preload3=os.path.join("salmon.dir", "salmon.transcripts.tpms.load")
    output:
        os.path.join("salmon.dir", "salmon.genes.tpms.load")
    params:
        job_mem=config["sqlite"].get("himem")
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("salmon.dir", "salmon.genes.tpms.load.log")
    run:
            statement = load(
                    infile=input.tpms,
                    outfile=str(output[0]),
                    job_memory={params.job_mem},
                    options="-i gene_id"
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

# --------------------- Run TXImeta --------------------------- #
rule tximeta:
    input:
        quant=expand(os.path.join("salmon.dir","{samp}","quant.genes.sf"), samp=S.samples.keys()),
    output:
        os.path.join("tximeta.dir","tximeta.RDS")
    params:
        geneset_path = os.path.join(config["txseq"].get("annotations"),
                                "api.dir/txseq.geneset.gtf.gz"),
        transcript_path = os.path.join(config["txseq"].get("annotations"),
                                "api.dir/txseq.transcript.fa.gz"),
        script=os.path.join(workflow.basedir,
                "R","scripts",
                "tximeta.R"),
        samples=config["samples"],
        tximeta_release=config["tximeta"].get("release"),
        tximeta_genome=config["tximeta"].get("genome"),
        tximeta_organism=config["tximeta"].get("organism"),
        txseq_salmon_index_name=config["txseq"].get("index_name")

    threads: 1
    resources:
        runtime=config["tximeta"].get("runtime"),
        mem_mb=config["tximeta"].get("memory"),
    log:
        os.path.join("tximeta.dir","tximeta.log")
    shell:
        """
        {{
            Rscript {params.script} \
                   --index_name={params.txseq_salmon_index_name} \
                   --salmondir=salmon.dir \
                   --transcripts={params.transcript_path} \
                   --geneset={params.geneset_path} \
                   --organism={params.tximeta_organism} \
                   --genomeversion={params.tximeta_genome} \
                   --release={params.tximeta_release} \
                   --samples={params.samples} \
                   --outfile={output} 
        }} > {log} 2>&1
        """

# --------------------- Get transcript info --------------------------- #
rule transcript_info:
    input:
        ti=os.path.join(config["txseq"].get("annotations"),
                          "api.dir/txseq.transcript.info.tsv.gz"),
        preload4=os.path.join("salmon.dir", "salmon.genes.tpms.load")
    output:
        os.path.join("transcript.info.load")
    params:
    threads: 1
    resources:
        **RESOURCES,
    log:
        "transcript.info.load.log"
    run:
            statement = load(
                    infile=input.ti,
                    outfile=str(output[0]),
                    options='-i "gene_id" -i "transcript_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

# --------------------- PGet n genes detected --------------------------- #
rule num_genes_detected:
    input:
       ti=os.path.join("transcript.info.load"),
       tpm=os.path.join("salmon.dir", "salmon.genes.tpms.load")
    output:
        os.path.join("qc.dir", "number.genes.detected.salmon.tsv.gz")
    params:
        database=config["database"].get("db"),
        script=os.path.join(workflow.basedir,
                "python",
                "salmon_no_genes_detected.py"),
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("qc.dir", "number.genes.detected.salmon.log")
    run:
            table = to_table(input.tpm)
            outfile=os.path.join("qc.dir", "number.genes.detected.salmon.tsv.gz")

            statement = f"""
                python {params.script} \
                   --database={params.database} \
                   --table={table} \
                   --outfile={outfile}
                """

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

## --------------------- Pload n genes into SQL --------------------------- #
rule load_n_genes_detected:
    input:
        ngd=os.path.join("qc.dir", "number.genes.detected.salmon.tsv.gz"),
        preload1=os.path.join("salmon.dir", "salmon.genes.load"),
        preload5=os.path.join("transcript.info.load")
    output:
        os.path.join("qc.dir","qc_no_genes_salmon.load")
    params:
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("qc.dir","qc_no_genes_salmon.load.log")
    run:
            statement = load(
                    infile=input.ngd,
                    outfile=str(output[0]),
                    options='-i "sample_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")