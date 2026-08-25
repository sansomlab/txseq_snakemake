# `Txseq`  ( Snakemake Implementation) 
> **⚠️ WARNING: `txseq`  is currently in pre-release and under active development. 
>Please use it with caution and report any issues by opening a GitHub issue.**
 
Txseq provides a suite of modular workflows for the QC and analysis of bulk RNA-sequencing data on high-performance 
compute clusters. These include pipelines for sanitising Ensembl annotations for RNA-seq analysis, building transcriptome 
indexes and performing quality control with Fastqc and Picard tools. Pipelines for traditional alignment-based mapping with 
HISAT2 and quantitation with `featureCounts` are provided along with a more-modern pseudo-alignment workflow that incorporates 
Salmon and tximeta. Template reports for sequence quality, post-mapping qc analysis, exploratory data analysis and differential 
gene expression are provided as Rmarkdown notebooks.


