# Mouse HSCs example

## Introduction

For this example we use paired end RNA-seq data from [Regan-Komito et al 'GM-CSF drives dysregulated hematopoietic stem cell activity and pathogenic extramedullary myelopoiesis in experimental spondyloarthritis' Nature Communicaitons 2020](https://doi.org/10.1038/s41467-019-13853-4). The data can be retrieved from the European Nucleotide Archive (ENA) under accession [PRJNA521342](https://www.ebi.ac.uk/ena/browser/view/PRJNA521342).

!!! Note 
    If you are working in the Kennedy workspace on the Oxford BMRC cluster, these data are already available in the "/well/kir/projects/mirror/ena/PRJNA521342/" folder.

## 1. Getting the configuration files and report templates

Make a suitable local folder and copy the samples.tsv and libraries.tsv for the "mouse_hscs" example into it. Here we choose to run the example in a folder with the path "~/work/hscs_example"

```text
  mkdir ~/work/txseq_hscs_example
  cd ~/work/txseq_hscs_example
  cp -r /path/to/txseq/examples/mouse_hscs/* .
  cp -r /path/to/txseq/reports .
```

Edit the txseq/libraries.tsv file to point to the location of the FASTQ files retrieved from the ENA on your system.

## 2. Retrieving and preparing annotations 

!!! Note 
    If you are working in the Kennedy workspace on the Oxford BMRC cluster, the required txseq annotations are available in the mirror folder and this step can be skipped.
  
Follow the instructions on the [Genomes and Annotations](genomesAndAnnotations.md) page to:

* Retrieve the required sequence and annotation files.

* Run the txseq ensembl pipeline to prepare santitised mouse genome and annotation files for RNA-seq analysis.

  
## 3 Building the Salmon and Hisat2 transcriptome indexes

!!! Note 
    If you are working in the Kennedy workspace on the Oxford BMRC cluster, txseq mouse transcriptome indexes have already been built in the mirror folder and this step can be skipped.

The example requires Salmon and Hisat2 indexes. To build these, follow the [Building transcriptome indexes](buildingTranscriptomeIndexes.md) instructions.


## 4. Checking FASTQ read quality

First we run the [FastQC tool](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) using the fastqc.smk pipeline. Get a copy of the configuration file and edit appropriately to ensure that the correct "libraries.tsv" and "samples.tsv" files paths are set.

```text
  module purge
  unset PYTHONPATH
  source ./txseq_smk_venv/bin/activate
  module load FastQC/0.12.1-Java-11
  cd ~/work/txseq_hscs_example
  mkdir fastqc
  cd fastqc
  txseq fastqc config # get a copy of the default configuration file
  emacs fastqc.yaml # edit the configuration file as appropriate
```

The pipeline can then be run as follows:

```text
  txseq fastqc full --jobs 20 
```

The pipeline parses and store the FastQC output into sqlite database in a file called "csvdb". To visualise the results you can use the associated R Markdown Report ("reports/fastqc.Rmd") in Rstudio.


## 5. Mapping with Hisat2

Next we map the data using the hisat.smk pipeline. Fetch and edit a copy of the configuration file to set the paths to the "libraries.tsv" and "samples.tsv" files and hisat index:

```text
  module purge
  unset PYTHONPATH
  source ./txseq_smk_venv/bin/activate
  module load HISAT2/2.2.1-gompi-2023a
  module load SAMtools/1.18-GCC-12.3.0
  cd ~/work/txseq_hscs_example
  mkdir hisat
  cd hisat
  txseq hisat config # get a copy of the default configuration file
  emacs pipeline_hisat.yml # edit the configuration file as appropriate
```
  
The pipeline can then be run as follows:
```text
  txseq hisat full --jobs 20
```
The output BAM files are located in the "hisat.dir" sub-directory.

## 6. Generating post-mapping QC statistics

After mapping with Hisat2, post-mapping QC statistics are computed using the bamqc.smk pipeline. This pipeline runs several [Picard](https://broadinstitute.github.io/picard/) tools including *CollectRnaSeqMetrics*, *EstimateLibraryComplexity*, *AlignmentSummaryMetrics* and *CollectInsertSizeMetrics* as well as some custom scripts. 

```text
  module purge
  module load R/4.5.1-gfbf-2023a-bare-noSciPy 
  unset PYTHONPATH
  source ./txseq_smk_venv/bin/activate
  module load picard/2.25.1-Java-11
  module load SAMtools/1.18-GCC-12.3.0
  module load Kent_tools/468-GCC-12.3.0
  cd ~/work/txseq_hscs_example
  mkdir bamqc
  cd bamqc
  txseq bamqc config # get a copy of the default configuration file
  emacs pipeline_bamqc.yml # edit the configuration file as appropriate
```

**IMPORTANT** You must unset the PYTHONPATH after loading the R module. 

The pipeline can then be run as follows:
```text
  txseq bamqc full --jobs 50
```
The results are saved in an sqlite database in the "csvdb" file. 

## 7. Quantitation with FeatureCounts

Count tables can be extracted from the BAM file using the feature_counts.smk pipeline. 

```text
  module purge
  unset PYTHONPATH
  source ./txseq_smk_venv/bin/activate
  module load Subread/2.0.6-GCC-12.3.0
  cd ~/work/txseq_hscs_example
  mkdir feature_counts
  cd feature_counts
  txseq feature_counts config # get a copy of the default configuration file
  emacs pipeline_feature_counts.yml # edit the configuration file as appropriate
```
The pipeline can then be run as follows:
```text
  txseq feature_counts full --jobs 20
```
The results are saved in an sqlite database in the "csvdb" file. 

## 8. Quantitation with Salmon

To quantitate the data using the salmon.smk pipeline, we begin by fetching and edit a copy of the configuration file to set the paths to the "libraries.tsv" and "samples.tsv" files and salmon index:
```text
  module purge
  module load R/4.5.1-gfbf-2023a-bare-noSciPy
  unset PYTHONPATH
  module load Salmon/1.11.4-GCC-12.3.0
  source ./txseq_smk_venv/bin/activate
  cd ~/work/txseq_hscs_example
  mkdir salmon
  cd salmon
  txseq salmon config # get a copy of the default configuration file
  emacs pipeline_salmon.yml # edit the configuration file as appropriate
```
**IMPORTANT** You must unset the PYTHONPATH after loading the R module. 
The pipeline can then be run as follows :
```text
  txseq salmon make full -v5 -p20
```
The results of the pipeline are stored in the "csvdb" sqlite database and as a tximeta object in the "tximeta.dir/tximeta.RDS" for downstream analysis. Flat tables of TPMs can be retrieved from the database or from the "salmon.dir/salmon.transcripts.tpms.txt.gz" file.

## 9. Post-mapping QC analysis

After running the bamqc and salmon pipelines post-mapping QC can be performed using the "post_mapping_qc.Rmd" report template.

Make a copy of the Rmd template file and open it in Rstudio to perform the analysis. The report visualises the individual QC statistics and performs a correlation analysis of the QC statistics with gene-expression space principle-components.

This analysis helps to identify confounding technical sources of variation.

## 10. Exploratory analysis

After running the salmon.smk pipeline the similarity between the samples in gene-expression space can be explored using the "exploratory_data_analysis.Rmd" R Markdown report template.

Make a copy of this file and open it in Rstudio to perform the analysis. The report produces plots showing hierarchical clustering of the samples by correlation of their expression profiles, the results of principle components analysis and a UMAP project of the samples.

Together with the post-mapping QC report this analysis is useful for the identification of outliers.

## 11. DESeq2 analysis

After running the salmon pipeline differential expression analysis can be performed using the "differential_expression.Rmd" R Markdown report template.

Make a copy of this file and open it in Rstudio to perform the analysis. 