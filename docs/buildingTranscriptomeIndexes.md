# Building transcriptome indexes

## Introduction


Txseq uses [Hisat2](http://daehwankimlab.github.io/hisat2/) for alignment-based gene expression quantification and [Salmon](https://github.com/COMBINE-lab/salmon) for quasi-alignment based gene expression quantification. To use these tools, it is first necessary to build method-specific transcriptomes indexes. 

As with preparation of the sanitised genome sequences and gene annotations, it is recommended to build transcriptome indexes in a central location for use in multiple projects.

Txseq has dedicated pipelines for building Salmon and Hisat2 indexes which can be used as follows.

!!! Note
    If you are using the KIR BMRC workspace, Salmon and Hisat index built with txseq can be found in the "/well/kir/projects/mirror/txseq/" directory.

## Building a Salmon Index

In a suitably named directory, load the required modules and obtain the salmon_index configuration file with the following commands:

```text
    module purge
	unset PYTHONPATH
	source ./txseq_smk_venv/bin/activate
	module load Salmon/1.11.4-GCC-12.3.0
	mkdir salmon_index
	cd salmon_index
	txseq salmon_index config
	emacs salmon_index.yaml 
```
After editing the yaml to point to the location of the "txseq.genome.fa.gz" and "txseq.transcript.fa.gz" files (see the "Genomes and Annotations" section), and configuring the parameters as desired, the pipeline can be executed with the following command:
```text
    txseq salmon_index full 
```

## Building a Hisat2 Index


In a suitably named directory, loads the required modules obtain the hisat_index configuration file with the following commands:
```text
    module purge
    unset PYTHONPATH
    source ./txseq_smk_venv/bin/activate
    module load HISAT2/2.2.1-gompi-2023a
    mkdir hisat_index
    cd hisat_index
    txseq hisat_index config
    emacs hisat_index.yaml
```
After editing the yaml to point to the location of the "txseq.genome.fa.gz" and "txseq.geneset.gtf.gz" files (see the "Genomes and Annotations" section), and configuring the parameters as desired, the pipeline can be executed with the following command:
```text
    txseq hisat_index full --jobs 2
```
