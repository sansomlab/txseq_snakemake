import os
from pathlib import Path
import re
import argparse
import logging
import sys
import gzip
import copy
import pandas as pd
import numpy as np

# <------------------------------ Logging ------------------------------------>

L = logging.getLogger(__name__)
log_handler = logging.StreamHandler(sys.stdout)
log_handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
log_handler.setLevel(logging.INFO)
L.addHandler(log_handler)
L.setLevel(logging.INFO)

# <------------------------------ Arguments ---------------------------------->

L.info("parsing arguments")

parser = argparse.ArgumentParser()
parser.add_argument("--outfile", default="None", type=str,
                    help='')
parser.add_argument("--coverage_histogram", default="None", type=str,
                    help='')


args = parser.parse_args()

L.info("Running with arguments:")
print(args)

# See comment above re histograms not always being generated
try:
        df = pd.read_csv(args.coverage_histogram, sep="\t")

        x = "normalized_position"
        cov = "All_Reads.normalized_coverage"

        three_prime_coverage = np.mean(df[cov][(df[x] > 70) & (df[x] < 90)])
        transcript_body_coverage = np.mean(df[cov][(df[x] > 20) & (df[x] < 90)])
        bias = three_prime_coverage / transcript_body_coverage

except:
        bias = np.nan

with open(args.outfile, "w") as out_file:
        out_file.write("three_prime_bias\n")
        out_file.write("%.2f\n" % bias)