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

rule full:
    input:
        genome_fasta=config["genome_fasta"],
        transcript_fasta=config["transcript_fasta"]
    output:
        os.path.join(config["rundir"],"salmon.index")
    params:
        kmerLens=config["kmerLen"],
        options=config["options"],
        sources=config["fasta_source"]
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
        }} > {log} 2>&1
        """
