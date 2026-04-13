import os
from pathlib import Path
import re
import argparse
import logging
import sys
import gzip
import copy

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
parser.add_argument("--pa", default=None, type=str,
                    help=(""))
parser.add_argument("--rundir", default="None", type=str,
                    help='')
parser.add_argument("--usedir", default="None", type=str,
                    help='')
parser.add_argument("--txfa", default="None", type=str,
                    help='')
parser.add_argument("--gtf", default=None, type=str,
                    help=(""))
parser.add_argument("--tx2gene", default=None, type=str,
                    help=(""))
parser.add_argument("--txinfo", default=None, type=str,
                    help=(""))
parser.add_argument("--outfile", default=None, type=str,
                    help=("name of the gzip compressed outfile"))

args = parser.parse_args()

L.info("Running with arguments:")
print(args)





Path(args.usedir).mkdir(parents=True, exist_ok=True)

                
print(args.pa)
 
if not os.path.exists(args.pa):
    raise ValueError("ypar masked primary assembly file not found")



src=os.path.join(args.rundir,args.pa)
print(src)
dest=os.path.join(args.usedir,"txseq.genome.fa.gz")
print(dest)

if os.path.exists(dest) or os.path.islink(dest):
    os.unlink(dest)

os.symlink(src,dest)
    

            
if not os.path.exists(args.txfa):
    raise ValueError("filtered transcript fasta file not found")

src=os.path.join(args.rundir,args.txfa)
print(src)
dest=os.path.join(args.usedir,"txseq.transcript.fa.gz")
print(dest)

if os.path.exists(dest) or os.path.islink(dest):
    os.unlink(dest)

os.symlink(src,dest)
                


if not os.path.exists(args.gtf):
    raise ValueError("filtered geneset gtf file not found")

src=os.path.join(args.rundir,args.gtf)
print(src)
dest=os.path.join(args.usedir,"txseq.geneset.gtf.gz")
print(dest)

if os.path.exists(dest) or os.path.islink(dest):
    os.unlink(dest)

os.symlink(src,dest)
     
            
if not os.path.exists(args.tx2gene):
   raise ValueError("transcript-to-gene map file not found")

src=os.path.join(args.rundir,args.tx2gene)
print(src)
dest=os.path.join(args.usedir,"txseq.transcript.to.gene.map")
print(dest)

if os.path.exists(dest) or os.path.islink(dest):
    os.unlink(dest)

os.symlink(src,dest)
                

            
if not os.path.exists(args.txinfo):
   raise ValueError("transcript info file not found")

src=os.path.join(args.rundir,args.txinfo)
print(src)
dest=os.path.join(args.usedir,"txseq.transcript.info.tsv.gz")
print(dest)

if os.path.exists(dest) or os.path.islink(dest):
    os.unlink(dest)

os.symlink(src,dest)
                
Path(args.outfile).touch()