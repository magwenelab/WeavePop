# Configuration

To configure the workflow you need to provide input files, edit the configuration file `config/config.yaml` and the execution profile `config/default/config.yaml` (for local execution) or `config/slurm/config.yaml` (for SLURM execution).

In the descriptions below when it mentions "specified in `field:`" it means a field in the `config/config.yaml`. When it is written as `field1: field2:` the second one is nested.

## Input files and their configuration

### FASTQ files:

Paired-end short-read FASTQ files, one forward and one reverse file for each sample. The names of these files should be the names used in the metadata `sample` column, followed by an extension specified in `fastq_suffix1:` and `fastq_suffix2:`. The files for all samples should be in one directory specified in `fastqs_directory:`. They can be gzip compressed.

### Reference genomes:

The names of the files must be the ones in the `lineage` column of the metadata (e.g. `VNI.fasta` and `VNI.gff`). They should be in one directory specified in `references_directory:`.

-   Providing annotations: Specify it with `annotate_references: activate: False`. Provide the FASTA and GFF files for each reference genome.

-   Annotating with main reference: Specify it with `annotate_references: activate: True`. If you want to have a common naming scheme for the genes or don't have an annotation (GFF file) of all your reference genomes you can provide one annotated main reference to annotate the rest using Liftoff. For this, provide the FASTA and GFF files for the main reference and specify the file names in `annotate_references: fasta:` and `annotate_references: gff:`. And provide only the FASTA file for each of the other reference genomes. If you activate the annotation of the reference genomes, all of them will be annotated.

### Metadata:

A comma-separated table with one sample per row. [Example](https://github.com/magwenelab/DiversityPipeline/blob/main/test/config/metadata.csv). Path specified in `metadata:`.\
Mandatory columns with these exact names:

-   `sample`: sample ID used in the FASTQ file names (no special characters or spaces).
-   `lineage`: lineage or group name that associates the sample with a reference genome (no special characters or spaces).
-   `strain`: strain name (a "common name" for each sample, it can be the same as `sample` if you don't have a different one).
-   If the plotting will be activated, you need one metadata column to color the samples. Specify the name of this column in `plotting: metadata2color:`. More columns with free format are allowed.

### Chromosomes:

A comma-separated table with one row per chromosome per lineage. [Example](https://github.com/magwenelab/DiversityPipeline/blob/main/test/config/chromosomes.csv). Path specified in `chromosomes:`.\
Mandatory columns with these exact names:

-   `lineage`: Lineage name (the same as in the metadata table and the names of the reference files).
-   `accession`: Sequence ID of the chromosomes in the FASTA and GFF of the reference of each lineage. Make sure each chromosome ID is not repeated in this file.
-   `chromosome`: Common name of the chromosome, e.g. chr01, 1, VNI_chr01.

### Exclude samples (optional)

If you want to exclude some of the samples in your metadata file from all analyses, you can provide a file with a list of sample names to exclude. Without a column name. Specify its path in `samples_to_exclude:`.

### Repeats database (optional)

Database of repetitive sequences in FASTA format to use for RepeatMasker. Needed for the CNV, plotting, and database modules. If you don't need a good identification of repeats specify it with `use_fake_database: True` and don't provide this file. We recommend the [RepBase database](https://www.girinst.org/server/RepBase/). You need to download it, extract the files, and concatenate them all in one FASTA file `config/RepBase_<version>.fasta`. Specify its path in `repeats_database:`.

```         
# Download the latest version and run the following:
tar -xvzf RepBase<version>.fasta.tar.gz
cat RepBase<version>.fasta/*.ref > RepBase.fasta
cat RepBase<version>.fasta/appendix/*.ref >> RepBase_<version>.fasta
rm -rf RepBase<version>.fasta/ RepBase<version>.fasta.tar.gz
```

### Genetic features for plots (optional)

If you want genetic features to be plotted in the depth and MAPQ plots, provide a comma-separated table with one row per gene. [Example](https://github.com/magwenelab/DiversityPipeline/blob/main/test/config/loci.csv). Specify its path in `plotting: loci:`.

Mandatory columns with these exact names:

-    `gene_id`: with the gene IDs (IDs of the reference genome's GFF).

-    `feature`: name of the feature (locus, pathway, centromere, individual gene name, etc.) the gene belongs to. Max 8 features.
