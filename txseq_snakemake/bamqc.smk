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

from txseq_snakemake.functions_from_cgat import pipeline as P

from txseq_snakemake.functions_from_cgat.pipeline.database import concatenate_and_load

from txseq_snakemake.functions_from_cgat.pipeline.parameters import PARAMS

import txseq_snakemake.functions_from_cgat.tables 
import txseq_snakemake.functions_from_cgat.csv2db 

PARAMS.update(config)

RESOURCES = config.get("resources", {})

## use txseq tasks to read in sample and library tables
S = T.samples(sample_tsv = config.get("samples"))

#get bam file name
def bam_file(sample_id):
    bam_fl=os.path.join(config.get("bam_path"),sample_id+".bam")
    return(bam_fl)

strand=config.get("strand")
print(strand)

paired=config.get("paired")
estimate_lib_comp=config["run"].get("estimateLibraryComplexity")

#get picard strand
if strand=="forward":
        strand_param = "FIRST_READ_TRANSCRIPTION_STRAND"
elif strand=="reverse":
        strand_param = "SECOND_READ_TRANSCRIPTION_STRAND"
else:
    strand_param = "NONE"
print(strand_param)

#get collectrnaseqmetrics options
if config["picard"].get("collectrnaseqmetrics_options"):
        picard_options_rnaseqmet = config["picard"].get("collectrnaseqmetrics_options")
else:
        picard_options_rnaseqmet = ""

#get alignment_summary_metrics options
if config["picard"].get("alignmentsummarymetric_options"):
        picard_options_alignsummet = config["picard"].get("alignmentsummarymetric_options")
else:
        picard_options_alignsummet = ""

#get estimate library complexity options
if config["picard"].get("estimatelibrarycomplexity_options"):
        picard_options_estlibcomp = config["picard"].get("estimatelibrarycomplexity_options")
else:
        picard_options_estlibcomp = ""

#get einsert size metrics options
if config["picard"].get("picard_insertsizemetric_options"):
        picard_options_insertsize = config["picard"].get("picard_insertsizemetric_options")
else:
        picard_options_insertsize = ""

#establish whether to run estimate library complexity

if paired and estimate_lib_comp:
    elc_torun=expand(os.path.join("bam.qc.dir/estimate.library.complexity.dir/","{samp}.library.complexity"),samp=S.samples.keys())
else:
    elc_torun=""

if paired and estimate_lib_comp:
    elc_torun_load=os.path.join("bam.qc.dir", "estimate_library_complexity.load")
else:
    elc_torun_load=""

#establish whether to run insert size metrics

if paired:
    insert_torun=[expand(os.path.join("bam.qc.dir/insert.size.metrics.dir/",
                                "{samp}.insert.size.metrics.summary"),samp=S.samples.keys()),
                    expand(os.path.join("bam.qc.dir/insert.size.metrics.dir/",
                                "{samp}.insert.size.metrics.histogram.pdf"),samp=S.samples.keys()),
                    expand(os.path.join("bam.qc.dir/insert.size.metrics.dir/",
                                "{samp}.insert.size.metrics.histogram"),samp=S.samples.keys())]
else:
    insert_torun=""



#rule all
rule full:
    input:
        rna_seq_metrics=expand(os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.rnaseq.metrics"), samp=S.samples.keys()),
        three_prime_bias=expand(os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.three.prime.bias"),samp=S.samples.keys()),
        align_sum_met=expand(os.path.join("bam.qc.dir", "alignment.summary.metrics.dir","{samp}.alignment.summary.metrics"),samp=S.samples.keys()),
        fs=expand(os.path.join("bam.qc.dir/fraction.spliced.dir/","{samp}.fraction.spliced"),samp=S.samples.keys()),
        elc=elc_torun,
        ins=insert_torun,
        rna_seq_metrics_load=os.path.join("bam.qc.dir", "qc_rnaseq_metrics.load"),
        three_prime_bias_load=os.path.join("bam.qc.dir", "three_prime_bias.load"),
        elc_load=elc_torun_load

## first pass of hisat. 
rule geneset_flat:
    input:
        os.path.join(config["txseq"].get("annotations"),
                            "api.dir/txseq.geneset.gtf.gz")
    output:
        os.path.join("annotations.dir","geneset.flat.gz")
    params:
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("annotations.dir","geneset.flat.log")
    shell:
        """
        {{
            gtfToGenePred \
                    -genePredExt \
                    -geneNameAsName2 \
                    -ignoreGroupsWithoutExons \
                    {input} \
                    /dev/stdout | \
                    awk 'BEGIN {{ OFS="\\t"}} \
                         {{print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}' \
                    | gzip -c \
                    > {output}
        }} > {log} 2>&1
        """

# ------------------- Picard: CollectRnaSeqMetrics -------------------------- #

rule collectRnaSeqMetrics:
    input:
        os.path.join("annotations.dir","geneset.flat.gz")
    output:
        os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.rnaseq.metrics")
    params:
        bam=lambda wildcards: bam_file(wildcards.samp),
        validation_stringency = config["picard"].get("validation_stringency"),
        coverage_out = os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.rnaseq.metrics.cov.hist"),
        chart_out = os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.rnaseq.metrics.cov.pdf"),
        mktemp_template = "ctmp.CollectRnaSeqMetrics.XXXXXXXXXX",
        picard_cmd=config["picard"].get("cmd"),
        strand_param=strand_param,
        picard_options=picard_options_rnaseqmet

    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.rnaseq.metrics.log")
    shell:
        """
        {{
            picard_out=`mktemp -p . {params.mktemp_template}`;
                   {params.picard_cmd} CollectRnaSeqMetrics \
                   -I {params.bam} \
                   --REF_FLAT {input} \
                   -O $picard_out \
                   --CHART {params.chart_out} \
                   --STRAND_SPECIFICITY {params.strand_param} \
                   --VALIDATION_STRINGENCY {params.validation_stringency} \
                   {params.picard_options};
                   grep . $picard_out | grep -v "#" | head -n2 \
                   > {output};
                   grep . $picard_out \
                   | grep -A 102 "## HISTOGRAM" \
                   | grep -v "##" \
                   > {params.coverage_out} || touch {params.coverage_out};
                   rm $picard_out;
        }} > {log} 2>&1
        """

# --------------------- load these RNAseq metrics --------------------------- #

rule collectRnaSeqMetrics_load:
    input:
        expand(os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.rnaseq.metrics"), samp=S.samples.keys()),
    output:
        os.path.join("bam.qc.dir", "qc_rnaseq_metrics.load")
    params:
    threads: 1
    resources:
        **RESOURCES,
    log:
        "bam.qc.dir/qc_rnaseq_metrics.log"
    run:
            statement = concatenate_and_load(
                    infiles=list(input),
                    outfile=str(output[0]),
                    regex_filename=r".*/.*/(.*)\.rnaseq\.metrics$",
                    cat="sample_id",
                    job_memory="8G",
                    options='-i "sample_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")



# --------------------- Three prime bias analysis --------------------------- #

rule three_prime_bias:
    input:
        os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.rnaseq.metrics.cov.hist")
    output: 
        os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.three.prime.bias")
    params:
        script=os.path.join(workflow.basedir,
                "python",
                "three_prime_bias.py"),
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.three.prime.bias.log")
    shell:
        """
        {{
            python {params.script} \
                --outfile={output} \
                --coverage_histogram={input}

        }} > {log} 2>&1
        """



# --------------------- load 3' bias --------------------------- #
rule three_prime_bias_load:
    input:
        expand(os.path.join("bam.qc.dir", "rnaseq.metrics.dir","{samp}.three.prime.bias"),samp=S.samples.keys()),
    output:
        os.path.join("bam.qc.dir", "three_prime_bias.load")
    params:
    threads: 1
    resources:
        **RESOURCES,
    log:
        "bam.qc.dir/three_prime_bias.log"
    run:
            statement = concatenate_and_load(
                    infiles=list(input),
                    outfile=str(output[0]),
                    regex_filename=".*/.*/(.*).three.prime.bias",
                    cat="sample_id",
                    job_memory="8G",
                    options='-i "sample_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

# --------------------- Alignment summary metrics --------------------------- #


rule alignment_summary_metrics:
    input:
        reference_sequence = os.path.join(config["txseq"].get("annotations"),
                                      "api.dir/txseq.genome.fa.gz")
    output:
        os.path.join("bam.qc.dir", "alignment.summary.metrics.dir","{samp}.alignment.summary.metrics")
    params:
        bam=lambda wildcards: bam_file(wildcards.samp),
        picard_options=picard_options_alignsummet,
        validation_stringency = config["picard"].get("validation_stringency"),
        mktemp_template = "ctmp.CollectAlignmentSummaryMetrics.XXXXXXXXXX",
        picard_cmd=config["picard"].get("cmd"),
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("bam.qc.dir", "alignment.summary.metrics.dir","{samp}.alignment.summary.metrics.log")
    shell:
        """
        {{
            picard_out=`mktemp -p . {params.mktemp_template}`;
                   {params.picard_cmd} CollectAlignmentSummaryMetrics \
                   -I {params.bam} \
                   -O $picard_out \
                   --REFERENCE_SEQUENCE {input} \
                   --VALIDATION_STRINGENCY {params.validation_stringency} \
                   {params.picard_options};
                   sed -e '1,/## HISTOGRAM/!d' $picard_out \
                   | grep . | grep -v "#" \
                   > {output};
                   rm $picard_out;
        }} > {log} 2>&1
        """


# --------------------- Fraction of spliced reads --------------------------- #


rule fraction_spliced:
    input:
        reference_sequence = os.path.join(config["txseq"].get("annotations"),
                                      "api.dir/txseq.genome.fa.gz")
    output:
        os.path.join("bam.qc.dir/fraction.spliced.dir/","{samp}.fraction.spliced")
    params:
        bam=lambda wildcards: bam_file(wildcards.samp),
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("bam.qc.dir/fraction.spliced.dir/","{samp}.fraction.spliced.log")
    shell:
        """
        {{
            echo "fraction_spliced" > {output};
                   samtools view {params.bam} \
                   | grep NH:i:1 \
                   | cut -f 6 \
                   | awk '{{if(index($1,"N")==0){{us+=1}} \
                           else{{s+=1}}}} \
                          END{{print s/(us+s)}}' \
                    >> {output}
        }} > {log} 2>&1
        """

# --------------------- Estimate library complexity --------------------------- #


rule estimate_library_complexity:
    input:
        reference_sequence = os.path.join(config["txseq"].get("annotations"),
                                      "api.dir/txseq.genome.fa.gz")
    output:
        os.path.join("bam.qc.dir/estimate.library.complexity.dir/",
                                "{samp}.library.complexity")
    params:
        bam=lambda wildcards: bam_file(wildcards.samp),
        mktemp_template = "ctmp.EstimateLibraryComplexity.XXXXXXXXXX",
        picard_cmd=config["picard"].get("cmd"),
        validation_stringency = config["picard"].get("validation_stringency"),
        picard_options=picard_options_estlibcomp
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("bam.qc.dir/estimate.library.complexity.dir/",
                                "{samp}.library.complexity.log")
    shell:
        """
        {{
            picard_out=`mktemp -p . {params.mktemp_template}`;
                   {params.picard_cmd} EstimateLibraryComplexity \
                   -I {params.bam} \
                   -O $picard_out \
                   --VALIDATION_STRINGENCY {params.validation_stringency} \
                   {params.picard_options};
                   grep . $picard_out | grep -v "#" | head -n2 \
                   > {output};
                   rm $picard_out;
        }} > {log} 2>&1
        """


# --------------------- load estimate library complexity --------------------------- #
rule estimate_library_complexity_load:
    input:
        elc_torun
    output:
        os.path.join("bam.qc.dir", "estimate_library_complexity.load")
    params:
    threads: 1
    resources:
        **RESOURCES,
    log:
        "bam.qc.dir/estimate_library_complexity.log"
    run:
            statement = concatenate_and_load(
                    infiles=list(input),
                    outfile=str(output[0]),
                    regex_filename=".*/.*/(.*).library.complexity",
                    cat="sample_id",
                    job_memory="8G",
                    options='-i "sample_id"'
                )

            print("Running statement:\n", statement, flush=True)
            shell(f"{statement} > {log} 2>&1")

# --------------------- Insert size metrics and histograms --------------------------- #


rule insert_size_metrics_and_hists:
    input:
        reference_sequence = os.path.join(config["txseq"].get("annotations"),
                                      "api.dir/txseq.genome.fa.gz")
    output:
        summary=os.path.join("bam.qc.dir/insert.size.metrics.dir/",
                                "{samp}.insert.size.metrics.summary"),
        hist_pdf=os.path.join("bam.qc.dir/insert.size.metrics.dir/",
                                "{samp}.insert.size.metrics.histogram.pdf"),
        hist=os.path.join("bam.qc.dir/insert.size.metrics.dir/",
                                "{samp}.insert.size.metrics.histogram")
    params:
        bam=lambda wildcards: bam_file(wildcards.samp),
        mktemp_template = "ctmp.CollectInsertSizeMetrics.XXXXXXXXXX",
        picard_cmd=config["picard"].get("cmd"),
        validation_stringency = config["picard"].get("validation_stringency"),
        picard_options=picard_options_insertsize
    threads: 1
    resources:
        **RESOURCES,
    log:
        os.path.join("bam.qc.dir/insert.size.metrics.dir/",
                                "{samp}.insert.size.metrics.log")
    shell:
        """
        {{
            picard_out=`mktemp -p . {params.mktemp_template}`;
                   {params.picard_cmd} CollectInsertSizeMetrics \
                   -I {params.bam} \
                   -O $picard_out \
                   --Histogram_FILE {output.hist_pdf} \
                   --VALIDATION_STRINGENCY {params.validation_stringency} \
                   --REFERENCE_SEQUENCE {input.reference_sequence} \
                   {params.picard_options};
                   grep "MEDIAN_INSERT_SIZE" -A 1 $picard_out \
                   > {output.summary};
                   sed -e '1,/## HISTOGRAM/d' $picard_out \
                   > {output.hist};
                   rm $picard_out;
        }} > {log} 2>&1
        """
