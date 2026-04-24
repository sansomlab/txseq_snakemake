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
RESOURCES = config.get("resources", {})

## use txseq tasks to read in sample and library tables
S = T.samples(sample_tsv = config.get("samples"),
                 library_tsv = config.get("libraries"))

def fastq_input(sample_id):
    sample = S.samples[sample_id]
    if sample.paired:
        fastq_input = "-1 " + ",".join(sample.fastq["read1"]) +\
                      " -2 " + ",".join(sample.fastq["read2"])
    else:
        fastq_input = "-r " + ",".join(sample.fastq["read1"])
    return(fastq_input)



salmon_quant_options = ''
if not config.get('salmon_quant_options') is None:
        salmon_quant_options = config.get('salmon_quant_options')
    
def libtype(sample_id):
    sample = S.samples[sample_id]
    return(sample.salmon_libtype)


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
        transcript_infoload=os.path.join("transcript.info.load")

rule quant:
    input:
    output:
        os.path.join("salmon.dir","{samp}","quant.sf")
    params:
        fastq_input=lambda wildcards: fastq_input(wildcards.samp),
        salmon_quant_options=salmon_quant_options,
        libtype=lambda wildcards: libtype(wildcards.samp),
        tx2gene = os.path.join(config["txseq"].get("annotations"),"api.dir/txseq.transcript.to.gene.map"),
        txseq_salmon_index=config["txseq"].get("salmon_index"),
        job_threads=config["salmon"].get("threads"),
        output_dir=os.path.join("salmon.dir","{samp}")
    threads: 1
    resources:
        **RESOURCES,
    log:
                os.path.join("salmon.dir","{samp}_quant.log")
    shell:
        """
        {{
            salmon quant -i {params.txseq_salmon_index} \
                                -p {params.job_threads} \
                                -g {params.tx2gene} \
                                {params.salmon_quant_options} \
                                -l {params.libtype} \
                                {params.fastq_input} \
                                -o {params.output_dir}
        }} > {log} 2>&1
        """

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
                    infiles=list(input),
                    outfile=str(output[0]),
                    regex_filename=".*/(.*)/quant.genes.sf",
                    cat="sample_id",
                    job_memory={params.job_mem},
                    options="-i Name -i sample_id"
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")


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
        **RESOURCES,
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


rule loadSalmonTranscriptTPMs:
    input:
       os.path.join("salmon.dir", "salmon.transcripts.tpms.txt.gz")
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
                    infile=input[0],
                    outfile=str(output[0]),
                    job_memory={params.job_mem},
                    options="-i transcript_id"
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

rule loadSalmonGeneTPMs:
    input:
       os.path.join("salmon.dir", "salmon.genes.tpms.txt.gz")
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
                    infile=input[0],
                    outfile=str(output[0]),
                    job_memory={params.job_mem},
                    options="-i gene_id"
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

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
        txseq_salmon_index=config["txseq"].get("salmon_index")

    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("tximeta.dir","tximeta.log")
    shell:
        """
        {{
            Rscript {params.script} \
                   --indexdir={params.txseq_salmon_index} \
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

rule transcript_info:
    input:
       os.path.join(config["txseq"].get("annotations"),
                          "api.dir/txseq.transcript.info.tsv.gz")
    output:
        os.path.join("transcript.info.load")
    params:
    threads: 1
    resources:
        **RESOURCES,
    log:
        "transcrit.info.load.log"
    run:
            statement = load(
                    infile=input[0],
                    outfile=str(output[0]),
                    job_memory="15G",
                    options='-i "gene_id" -i "transcript_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")