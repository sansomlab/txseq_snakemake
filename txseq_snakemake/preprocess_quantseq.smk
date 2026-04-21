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

fqs=pd.read_csv(config.get("libraries"),sep="\t")

SAMPLE_TO_FASTQ = dict(zip(fqs["sample_id"], fqs["fastq_path"]))
SAMPLES = fqs["sample_id"].tolist()


rule full:
    input:
        expand("cutadapt.dir/{sample}.trimmed.fq.gz", sample=SAMPLES)

rule cut:
    input:
        lambda wc: SAMPLE_TO_FASTQ[wc.sample]
    output:
        "cutadapt.dir/{sample}.trimmed.fq.gz"
    params:
        trim=config["trim"],
        mktemp_template = "ctmp.cutadapt.XXXXXXXXXX"
    threads: 1
    resources:
        **RESOURCES,
    log:
        "cutadapt.dir/{sample}.trimmed.log"
    shell:
        """
        {{
        if [ "{params.trim}" != "False" ] && [ "{params.trim}" != "false" ]; then
            echo "hi"
            trim_out=`mktemp -p . {params.mktemp_template}`;
            cutadapt -l {params.trim} {input} -o $trim_out; 
            cutadapt_out_a=`mktemp -p . {params.mktemp_template}`;
            cutadapt_out_b=`mktemp -p . {params.mktemp_template}`;
            cutadapt -m 20 -O 20 -a "polyA=A{{20}}" -a "QUALITY=G{{20}}" -n 2 $trim_out -o $cutadapt_out_a;
            cutadapt -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=A{{18}}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" $cutadapt_out_a -o $cutadapt_out_b;
            cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed $cutadapt_out_b -o {output};
            rm $cutadapt_out_a $cutadapt_out_b $trim_out;
        else
            echo "bye"
            trim_out={input};
            cutadapt_out_a=`mktemp -p . {params.mktemp_template}`;
            cutadapt_out_b=`mktemp -p . {params.mktemp_template}`;
            cutadapt -m 20 -O 20 -a "polyA=A{{20}}" -a "QUALITY=G{{20}}" -n 2 $trim_out -o $cutadapt_out_a;
            cutadapt -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=A{{18}}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" $cutadapt_out_a -o $cutadapt_out_b;
            cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed $cutadapt_out_b -o {output};
            rm $cutadapt_out_a $cutadapt_out_b;
        fi
        }} > {log} 2>&1
        """