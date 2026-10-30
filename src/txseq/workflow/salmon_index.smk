"""========================
salmon_index.smk
========================

Overview
--------

This pipeline makes a Salmon SAF genome index. This type of index uses the full genome as a decoy: according to the Salmon authors this type of index 'does the best job in avoiding spurious alignments to annotated transcripts'.

Configuration
-------------

The pipeline requires a configured :file:`salmon_index.yaml` file.


Input files
-----------

For building the SAF genome index the pipeline requires:

#. A gzip compressed FASTA file containing the transcript sequences.
#. A gzip compressed FASTA file containing the genome primary assembly sequences.

The location of these three files must be specified in the 'salmon_index.yaml' file.


Requirements
------------

The following software is required:

#. Salmon


Output files
------------

The pipeline generates a salmon index called "salmon_index" in the current folder.


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

#general pipeline prep
# Per-rule job resources (mem_mb, runtime, ...) from the `resources:` section
# of the config. "threads" is a directive, not a resource, so it is held
# separately. The Slurm account and partition are cluster-level settings
# that the txseq CLI applies to every job via --default-resources.
RESOURCES = dict(config.get("resources") or {})
THREADS = RESOURCES.pop("threads", 1)
source = "--gencode" if config["fasta_source"] == "gencode" else ""

# --------------------- Pipeline --------------------------- #
rule full:
    input:
        genome_fasta=config["genome_fasta"],
        transcript_fasta=config["transcript_fasta"]
    output:
        "salmon_index.sentinel"
    params:
        kmerLens=config["kmerLen"],
        options=config["options"],
        sources=source
        #salmon_index_loc=os.path.join(config["rundir"],"salmon.index")
    threads: 12
    resources:
        **RESOURCES,
    log:
        "index.log"
    shell:
        """
        {{
            grep "^>" <(gunzip -c {input.genome_fasta}) \
                   | cut -d " " -f 1 \
                   > decoys.txt;
                   sed -i.bak -e 's/>//g' decoys.txt;
                   cat {input.transcript_fasta} {input.genome_fasta} \
                   > gentrome.fa.gz;
                   salmon index -t gentrome.fa.gz \
                                -d decoys.txt \
                                -k {params.kmerLens} \
                                -p {threads} \
                                -i salmon_index \
                                {params.options} \
                                {params.sources}
                touch {output}
        }} > {log} 2>&1
        """
