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
        txinfo=os.path.join("ensembl.dir","api_completed.txt")

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
        txinfo=os.path.join("ensembl.dir","api_completed.txt")
    log:
        txinfo=os.path.join("ensembl.dir","api_completed.log")
    params:
    run:
            import os
            from pathlib import Path
            if not os.path.exists("ensembl.dir/api.dir"):
                os.mkdir("ensembl.dir/api.dir")

            pa = input.pa
    
            if not os.path.exists(pa):
                raise ValueError("ypar masked primary assembly file not found")
    
            os.symlink(os.path.join("..","..",pa), 
                    os.path.join("ensembl.dir","api.dir","txseq.genome.fa.gz"))
    
            txfa = input.txfa
            
            if not os.path.exists(txfa):
                raise ValueError("filtered transcript fasta file not found")
                
            os.symlink(os.path.join("..","..",txfa), 
                    os.path.join("ensembl.dir","api.dir","txseq.transcript.fa.gz"))
            
            gtf = input.gtf
            
            if not os.path.exists(gtf):
                raise ValueError("filtered geneset gtf file not found")
                
            os.symlink(os.path.join("..","..", gtf), 
                    os.path.join("ensembl.dir","api.dir","txseq.geneset.gtf.gz"))
            
            tx2gene = input.tx2gene
            
            if not os.path.exists(tx2gene):
                raise ValueError("transcript-to-gene map file not found")
                
            os.symlink(os.path.join("..","..", tx2gene), 
                    os.path.join("ensembl.dir","api.dir","txseq.transcript.to.gene.map"))
            
            txinfo = input.txinfo
            
            if not os.path.exists(txinfo):
                raise ValueError("transcript info file not found")
                
            os.symlink(os.path.join("..","..", txinfo), 
                    os.path.join("ensembl.dir","api.dir","txseq.transcript.info.tsv.gz"))
            
            touch output

