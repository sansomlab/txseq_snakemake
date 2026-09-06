"""===================
pipeline_ensembl.py
===================

Overview
--------

This pipeline post-processes Ensembl annotation files to prepare a set of reference files suitable for the analysis of RNA-sequencing data. We use here Ensembl sequences rather than Gencode sequences to avoid this issue: `Salmon issue 214 <https://github.com/COMBINE-lab/salmon/issues/214>`_ . Gencode and Ensembl GTF files are in any case equivalent.

This pipeline prepares genome and annotation files that (a) only include records from PRIMARY contigs, i.e. chromosomes, mitochondrial sequences, unlocalised and unplaced scaffolds but not alternative (ALT) contigs and (b) exclude multi-placed transcript sequences by masking/excluding sequence/records from the Y chromosome PAR region. 

It performs the following tasks:

#. Makes a version of primary assembly FASTA file in which the Y chromosome PAR regions are hard masked.
#. Makes a FASTA file containing the coding (cDNA) and non-coding (ncRNA) Ensembl transcripts from PRIMARY contigs that are not multi-placed transcript sequences.
#. Makes a similarly filtered version of the Ensembl geneset GTF file.
#. Makes a transcript -> gene map for use with sample
#. Makes a transcript information table that contains information on transcript and gene names and biotypes.


Configuration
-------------

The pipeline requires a configured :file:`pipeline_ensembl.yml` file.

Input files
-----------

The pipeline requires the following inputs

#. The Ensembl primary assembly FASTA sequences
#. The Ensembl geneset in GTF format
#. The Ensembl cDNA FASTA sequences
#. The Ensembl ncRNA FASTA sequences 
#. PAR region definitions in BED format

The location of these three files must be specified in the 'pipeline_ensembl.yml' file.

Output files
------------

The pipeline creates an "api" folder with the following files for use by downstream pipelines:

#. api.dir/txseq.genome.fa.gz: a copy of the Ensembl primary assembly in which Y PAR regions are hard masked
#. api.dir/txseq.transcript.fa.gz: all records from the Ensembl cDNA and ncRNA transcript fasta files that are on primary contigs and not in the Y PAR region
#. api.dir/txseq.geneset.gtf.gz: all records from the Ensembl gtf file that are on primary contigs and not in the Y PAR region
#. api.dir/txseq.transcript.to.gene.map: a tab-seperated list of transcript_id -> gene_id mappings for use with Salmon
#. api.dir/txseq.transcript.info.tsv.gz: a tab-seperated table of transcription information (including transcript_name, transcript_biotype, gene_name and gene_biotype)

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
import txseq.tasks as T

RESOURCES = config.get("resources", {})

# --------------------- Pipeline --------------------------- #
localrules: full
#rule all
rule full:
    input:
        os.path.join("ensembl.dir","api.dir","api.completed.txt")

# --------------------- Extract and hard mask chromosome Y PAR regions --------------------------- #

rule extract_YPAR:
    input: 
        config["par"]
    output: 
        os.path.join("ensembl.dir","Y.PAR.bed")
    resources:
        **RESOURCES,
    log: 
        os.path.join("ensembl.dir","Y.PAR.extract.log")
    shell:
        """
        {{
        grep Y {input} > {output}
        }} > {log} 2>&1
        """

rule hard_mask_YPAR:
    input:
        os.path.join("ensembl.dir","Y.PAR.bed")
    output:
        os.path.join("ensembl.dir","ypar.masked.primary.assembly.fa.gz")
    log: 
        os.path.join("ensembl.dir","Y.PAR.hard_mask.log")
    resources:
        **RESOURCES,
    params:
        primary_assembly=config["primary_assembly"],
        ypar_bed=os.path.join("ensembl.dir","Y.PAR.bed"),
        out=os.path.join("ensembl.dir","ypar.masked.primary.assembly.fa"),
    shell:
        """
        {{
        bedtools maskfasta \
                   -fi <(zcat {params.primary_assembly}) \
                   -fo {params.out} \
                   -bed {params.ypar_bed};
                gzip {params.out};
        }} > {log} 2>&1
        """

rule contigs:
    input: 
        os.path.join("ensembl.dir","ypar.masked.primary.assembly.fa.gz")
    output:
        os.path.join("ensembl.dir","contigs.txt")
    log: 
        os.path.join("ensembl.dir","contigs.log")
    resources:
        **RESOURCES,
    shell:
        """
        {{
            zgrep \> {input} \
            | sed 's/>//g' \
            > {output}
        }} > {log} 2>&1
        """


# --------------------- Filter transcripts --------------------------- #

rule filtered_transcript_fasta:
    input:
        contigs=os.path.join("ensembl.dir","contigs.txt"),
        ypar=os.path.join("ensembl.dir","Y.PAR.bed")
    output:
        os.path.join("ensembl.dir","filtered.transcripts.fa.gz")
    log:
        os.path.join("ensembl.dir","filtered.transcripts.log")
    resources:
        **RESOURCES,
    params:
        cdna=config["cdna"],
        ncrna=config["ncrna"],
        script=os.path.join(workflow.basedir,
                        "..",
                "scripts",
                        "ensembl_filter_transcript_fasta.py"),
    shell:
        """
        {{
            python {params.script} \
                 --ensembltxfasta={params.cdna},{params.ncrna} \
                 --contigs={input.contigs} \
                 --mask={input.ypar} \
                 --outfile={output}
        }} > {log} 2>&1
        """

# --------------------- Make filtered GTF --------------------------- #

rule filtered_gtf:
        input:
            contigs=os.path.join("ensembl.dir","contigs.txt"),
            ypar=os.path.join("ensembl.dir","Y.PAR.bed")
        output:
            os.path.join("ensembl.dir","filtered.geneset.gtf.gz")
        log:
            os.path.join("ensembl.dir","filtered.geneset.gtf.log")
        resources:
            **RESOURCES,
        params:
            script=os.path.join(workflow.basedir,
                        "..",
                "scripts",
                        "ensembl_filter_gtf.py"),
            geneset=config["geneset"],
        shell:
            """
            {{
                python {params.script} \
                    --ensemblgtf={params.geneset} \
                    --contigs={input.contigs} \
                    --mask={input.ypar} \
                    --outfile={output}
            }} > {log} 2>&1
            """

# --------------------- Make transcript to gene map --------------------------- #

rule transcript_to_gene_map:
    input: 
        os.path.join("ensembl.dir","filtered.geneset.gtf.gz")
    output:
        os.path.join("ensembl.dir","transcript_to_gene_map")
    log:
        os.path.join("ensembl.dir","transcript_to_gene_map.log")
    resources:
        **RESOURCES,
    shell:
        """
        {{
            zgrep transcript {input} \
            | sed 's/.*gene_id "\\([^"]*\\)".*transcript_id "\\([^"]*\\)".*$/\\2\\\t\\1/g' \
            | sort -u > {output}
        }} > {log} 2>&1
        """

# --------------------- Get transcript info --------------------------- #
rule transcript_info:
    input: 
        os.path.join("ensembl.dir","filtered.geneset.gtf.gz")
    output:
        os.path.join("ensembl.dir","transcript.info.tsv.gz")
    log:
        os.path.join("ensembl.dir","transcript.info.log")
    resources:
        **RESOURCES,
    params:
            script=os.path.join(workflow.basedir,
                        "..",
                "scripts",
                        "ensembl_extract_gtf_attributes.py"),
            geneset=config["geneset"],
    shell:
        """
        {{
            python {params.script} \
            --ensemblgtf={params.geneset} \
            --attributes=transcript_id,transcript_name,transcript_biotype,gene_id,gene_name,gene_biotype \
            --outfile={output}
        }} > {log} 2>&1
        """

# --------------------- Tidy up results --------------------------- #
rule api:
    input:
        pa=os.path.join("ensembl.dir","ypar.masked.primary.assembly.fa.gz"),
        txfa=os.path.join("ensembl.dir","filtered.transcripts.fa.gz"),
        gtf=os.path.join("ensembl.dir","filtered.geneset.gtf.gz"),
        tx2gene=os.path.join("ensembl.dir","transcript_to_gene_map"),
        txinfo=os.path.join("ensembl.dir","transcript.info.tsv.gz")
    output:
        os.path.join("ensembl.dir","api.dir","api.completed.txt")
    log:
        txinfo=os.path.join("ensembl.dir","api.completed.log")
    resources:
        **RESOURCES,
    params:
        rundir=config["rundir"],
        usedir=os.path.join(config["rundir"],"ensembl.dir","api.dir"),
        script=os.path.join(workflow.basedir,
                        "..",
                "scripts",
                        "ensembl_api.py"),
    shell:
            """
            {{
                python {params.script} \
                --pa={input.pa} \
                --rundir={params.rundir} \
                --txfa={input.txfa} \
                --gtf={input.gtf} \
                --tx2gene={input.tx2gene} \
                --txinfo={input.txinfo} \
                --usedir={params.usedir} \
                --outfile={output}
            }} > {log} 2>&1
            """

            