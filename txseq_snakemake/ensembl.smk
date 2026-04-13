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

#rule all
rule full:
    input:
        os.path.join("ensembl.dir","api_dir","api_completed.txt")

rule extract_YPAR:
    input: 
        config["par"]
    output: 
        os.path.join("ensembl.dir","Y.PAR.bed")
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
    shell:
        """
        {{
            zgrep \> {input} \
            | sed 's/>//g' \
            > {output}
        }} > {log} 2>&1
        """

rule filtered_transcript_fasta:
    input:
        contigs=os.path.join("ensembl.dir","contigs.txt"),
        ypar=os.path.join("ensembl.dir","Y.PAR.bed")
    output:
        os.path.join("ensembl.dir","filtered.transcripts.fa.gz")
    log:
        os.path.join("ensembl.dir","filtered.transcripts.log")
    params:
        cdna=config["cdna"],
        ncrna=config["ncrna"],
        script=os.path.join(workflow.basedir,
                        "python",
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

rule filtered_gtf:
        input:
            contigs=os.path.join("ensembl.dir","contigs.txt"),
            ypar=os.path.join("ensembl.dir","Y.PAR.bed")
        output:
            os.path.join("ensembl.dir","filtered.geneset.gtf.gz")
        log:
            os.path.join("ensembl.dir","filtered.geneset.gtf.log")
        params:
            script=os.path.join(workflow.basedir,
                        "python",
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

rule transcript_to_gene_map:
    input: 
        os.path.join("ensembl.dir","filtered.geneset.gtf.gz")
    output:
        os.path.join("ensembl.dir","transcript_to_gene_map")
    log:
        os.path.join("ensembl.dir","transcript_to_gene_map.log")
    shell:
        """
        {{
            zgrep transcript {input} \
            | sed 's/.*gene_id "\\([^"]*\\)".*transcript_id "\\([^"]*\\)".*$/\\2\\\t\\1/g' \
            | sort -u > {output}
        }} > {log} 2>&1
        """

rule transcript_info:
    input: 
        os.path.join("ensembl.dir","filtered.geneset.gtf.gz")
    output:
        os.path.join("ensembl.dir","transcript.info.tsv.gz")
    log:
        os.path.join("ensembl.dir","transcript.info.log")
    params:
            script=os.path.join(workflow.basedir,
                        "python",
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

rule api:
    input:
        pa=os.path.join("ensembl.dir","ypar.masked.primary.assembly.fa.gz"),
        txfa=os.path.join("ensembl.dir","filtered.transcripts.fa.gz"),
        gtf=os.path.join("ensembl.dir","filtered.geneset.gtf.gz"),
        tx2gene=os.path.join("ensembl.dir","transcript_to_gene_map"),
        txinfo=os.path.join("ensembl.dir","transcript.info.tsv.gz")
    output:
        os.path.join("ensembl.dir","api_dir","api_completed.txt")
    log:
        txinfo=os.path.join("ensembl.dir","api_completed.log")
    params:
        rundir=config["rundir"],
        usedir=os.path.join(config["rundir"],"ensembl.dir","api_dir"),
        script=os.path.join(workflow.basedir,
                        "python",
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

            