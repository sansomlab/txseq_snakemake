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

## use txseq tasks to read in sample and library tables
S = T.samples(sample_tsv = config.get("samples"))

rule full:
    input:
        expand(os.path.join("fastqc.dir", "{samp}"), samp=S.fastqs.keys())

rule fastqc:
    input:
        os.path.join(config["fq_dir"],"{samp}_1.fq.gz")
    output: 

