# Genomes and Annotations

## Introduction


In general, for RNA-seq analysis it is best to use genome sequences and gene annotations that:

* Include sequences and gene annotations from PRIMARY contigs, i.e. chromosomes, mitochondrial sequences, unlocalised and unplaced scaffolds

* Do not include sequences and gene annotations from alternative (ALT) contigs

* Exclude multi-placed transcript sequences by masking/excluding sequence/records from the Y chromosome PAR region. 

The txseq repository uses genome sequences and gene annotations retrieved from [Ensembl](https://www.ensembl.org/index.html). These are pre-processed using the "ensembl.smk" pipeline which addresses all the issues above and outputs sanitised genome sequence and annotation files.


!!! note-square "Note"
    If you are using the KIR BMRC workspace, sanitised genome sequences and gene annotations can be found in the <kbd>/well/kir/projects/mirror/txseq/</kbd> directory.



## Retrieving genome sequences and gene annotations

The following files are required:

* The Ensembl primary assembly FASTA sequences
* The Ensembl geneset in GTF format
* The Ensembl cDNA FASTA sequences
* The Ensembl ncRNA FASTA sequences 
* PAR region definitions in BED format


The current Ensembl genome and annotation files can retrieve from the [Ensembl FTP website](http://www.ensembl.org/info/data/ftp/index.html).

PAR region locations can be retrieved from e.g. the [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc/)


### Example 1: obtaining genome sequences and genes annotation files for analysis of human data

The human Ensembl genome and annotation files (for Ensembl release 115) can be retrieved using the following commands:

<div class="nord" markdown="1">
```py
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
```
</div>
    
The [human PAR coordinates](https://www.ncbi.nlm.nih.gov/grc/human) should then be used to prepare a bed file, for example for the GRCh38.p14 release of the human genome the (tab-separated) file should look like this:

    X	10001	2781479	PAR.1
    X	155701383	156030895	PAR.2
    Y	10001	2781479	PAR.1
    Y	56887903	57217415	PAR.2


### Example 2: Obtaining genome sequences and genes annotation files for analysis of mouse data

The mouse Ensembl genome and annotation files (for Ensembl release 115) can be retrieved using the following commands:

<div class="nord" markdown="1">
```py
wget https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz
```
</div>

The [mouse PAR locations](https://www.ncbi.nlm.nih.gov/grc/mouse) should thenbe used to prepare a bed file, for example for the GRCm39 release of the mouse genome the (tab-separated) file should look like this:

    X	168752755	169376592	PAR
    Y	90757114	91355967	PAR

## Preparing the txseq-sanitised genome and annotations

It is recommended to prepare txseq-sanitised genome sequences and annotations in a central location for use in all of your RNA-seq projects.

In a suitable directory, load the required modules and obtain a copy of the ensembl.smk configuration file:

<div class="nord" markdown="1">
```py
module purge
unset PYTHONPATH
source ./txseq_smk_venv/bin/activate
module load BEDTools/2.31.0-GCC-12.3.0
mkdir ensembl
cd ensembl
txseq ensembl config
emacs ensembl.yaml
```

    
After editing the .yml file to provide the locations of the Ensembl genome, Ensembl annotations and the PAR bed file, execute the pipeline with the following command:

```py
    txseq ensembl full
```
</div>
    
The output of the pipeline is an "api.dir" folder that contains the following files that can be used to build indexes for RNA-seq mapping and quantification tools:

* txseq.geneset.gtf.gz - the sanitised geneset
* txseq.genome.fa.gz - the sanitised and PAR masked genome
* txseq.transcript.fa.gz - the sanitised transcripts
* txseq.transcript.info.tsv.gz - a flat tsv table of transcript information (for the santitised transcript set)
* txseq.transcript.to.gene.map - a flat 2-column tsv table containing a map of transcripts -> genes (for the sanitised transcript set)