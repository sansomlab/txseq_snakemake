import sys
import shutil
import os
from pathlib import Path
import glob
import sqlite3

import sqlalchemy
from sqlalchemy import text

import pandas as pd
import numpy as np

# import local pipeline utility functions
import txseq_snakemake.tasks as T
import txseq_snakemake.tasks.readqc as readqc

from txseq_snakemake.functions_from_cgat.iotools import snip
import txseq_snakemake.functions_from_cgat.experiment as E
import txseq_snakemake.functions_from_cgat.iotools as IOTools
from txseq_snakemake.functions_from_cgat.pipeline.database import load 
from txseq_snakemake.functions_from_cgat.pipeline.parameters import PARAMS

PARAMS.update(config)
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


fastqc_outputs=[
    "adapter_content",
    "basic_statistics",
    "overrepresented_sequences",
    "per_base_n_content",
    "per_base_sequence_content",
    "per_base_sequence_quality",
    "per_sequence_gc_content",
    "per_sequence_quality_scores",
    "per_tile_sequence_quality",
    "sequence_duplication_levels",
    "sequence_length_distribution"
]

rule full:
    input:
        run_fastqc=expand(os.path.join("fastqc.dir", "{samp}.sentinel"), samp=S.fastqs.keys()),
        summarise_fastqc=expand("fastqc.summary.dir/fastqc_{module}.tsv.gz",module=fastqc_outputs),
        build_sum_status="fastqc.summary.dir/fastqc_status_summary.tsv.gz",
        all_stats_load=expand("fastqc.summary.dir/fastqc_{module}.load",
               module=fastqc_outputs),
        status_summary_stats_load="fastqc.summary.dir/fastqc_status_summary.load",
        load_meta="fastqc.dir/load.metadata.sentinel",

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


rule summariseFastQC:
    input:
        expand(os.path.join("fastqc.dir", "{samp}.sentinel"), samp=S.fastqs.keys())
    output: 
        expand("fastqc.summary.dir/fastqc_{module}.tsv.gz",
               module=fastqc_outputs)
    params:

    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "fastqc.summary.dir/summarise_fastqc.log"
    run:
        all_files = []

        for infile in input:

            track = snip(infile, ".sentinel")
            all_files.extend(glob.glob(
                os.path.join(track, "*_fastqc",
                            "fastqc_data.txt")))

        dfs = readqc.read_fastqc(
            all_files)
        
        print("Expected:", fastqc_outputs)
        print("Observed:", list(dfs.keys()))

        for key, df in dfs.items():
            fn = os.path.join("fastqc.summary.dir", f"fastqc_{key}.tsv.gz")
            E.info("writing to {}".format(fn))
            with IOTools.open_file(fn, "w") as outf:
                df.to_csv(outf, sep="\t", index=True)


rule buildFastQCSummaryStatus:
    input:
        expand(os.path.join("fastqc.dir", "{samp}.sentinel"), samp=S.fastqs.keys())
    output: 
        "fastqc.summary.dir/fastqc_status_summary.tsv.gz"
    params:
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "fastqc.summary.dir/fastqc_status_summary.log"
    run:
        readqc.buildFastQCSummaryStatus(
            list(input),
            output[0],
            "fastqc.dir")

rule loadFastQC:
    input:
        all_stats=expand("fastqc.summary.dir/fastqc_{module}.tsv.gz",
               module=fastqc_outputs),
        status_summary_stats="fastqc.summary.dir/fastqc_status_summary.tsv.gz",
    output: 
        all_stats_load=expand("fastqc.summary.dir/fastqc_{module}.load",
               module=fastqc_outputs),
        status_summary_stats_load="fastqc.summary.dir/fastqc_status_summary.load",
    params:
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "fastqc.summary.dir/load_fastqc.log",
    run:
        infiles=list(input.all_stats)+[input.status_summary_stats]
        outfiles = list(output.all_stats_load) + [output.status_summary_stats_load]
        for infile, outfile in zip(infiles, outfiles):
            print(infile)
            base = os.path.basename(infile).replace(".tsv.gz", "")
            print(outfile)
            print("here we goooooo")
        # a check to make sure file isn't empty
            n = 0
            with IOTools.open_file(infile) as f:
                for i, line in enumerate(f):
                    n =+ i
            if n > 0:
                print("this file")
                print(infile)
                print("has n>0")
                statement=load(str(infile), str(outfile), options="--add-index=track")
                print("Running statement:\n", statement, flush=True)
                shell(f"{statement} > {log} 2>&1")
            else:
                print("this file")
                print(infile)
                print("has n=0")
                table_name = base
                database_sql = config["database"].get("url")
                database_name = os.path.basename(database_sql)
                statement = f"""
                            sqlite3 {database_name} \
                            'DROP TABLE IF EXISTS {table_name}; \
                            CREATE TABLE {table_name} \
                            ("track" text PRIMARY KEY, "Sequence" text, \
                            "Count" integer, "Percentage" integer, "Possible Source" text);' \
                            'INSERT INTO {table_name} VALUES ("NA", "NA", 0, 0, "NA");'"""
                

                print("Running statement:\n", statement, flush=True)
                shell(f"{statement} > {log} 2>&1")




rule load_metadata:
    input:
    output: 
        "fastqc.dir/load.metadata.sentinel"
    params:
        sqlite_url=config["database"].get("url")
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    log:
        "fastqc.dir/load.metadata.log"
    run:
        
        db =sqlalchemy.create_engine(params.sqlite_url)

        with db.connect() as dbconn:

            S.sample_table.to_sql(name = 'samples',
                                con= dbconn, 
                                index= False, 
                                if_exists='replace') 
        
        

            dbconn.execute(text("CREATE INDEX samples_sample_id on samples (sample_id)"))
            
        with db.connect() as dbconn:
        
            S.fastq_table.to_sql(name = 'fastqs',
                                con= dbconn, 
                                index=False, 
                                if_exists='replace') 

            dbconn.execute(text("CREATE INDEX fastqs_fastq_id on fastqs (sample_id)"))
        
        db.dispose()

        IOTools.touch_file(str(output))

