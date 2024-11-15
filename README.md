
# Diversity Pipeline

## Description

This is a Snakemake workflow to map short-reads of samples to the reference genome of the corresponding lineage. You will get the mapping file, variant calling file, a reference-based assembly, and an annotation GFF. An SQL database will be created with the DNA and protein sequences, and with tables describing the presence of the called variants and their effects. Additionaly you can get analyses of the coverage and the mapping quality, plots with this analyses' results, including detection of copy-number variants, and a table with the intersection of the detected genetic variants between the samples of the same lineage.  

If you want to have a common naming scheme of your genes and/or you don't have GFF files of your reference genomes you can provide a main reference to lift over the annotation from this one to your reference genomes.

## Installation

### Install Conda

Install Mamba or Miniconda following the instructions from their webpage:
* Mamba: [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge) (recommended)
* Miniconda: [https://docs.anaconda.com/miniconda/](https://docs.anaconda.com/miniconda/)

After successfully installing conda, set strict channel priority by running  
`conda config --set channel_priority strict`

### Download the workflow 

#### Option 1: Download this GitHub repository

Use the green button `<> Code` and click `Download ZIP`.  
Extract the content of the `.zip` file.

#### Option 2: Use Snakedeploy (PENDING)

#### Option 3: Download from Bioconda (PENDING)


### Install the Snakemake Conda environment

In your terminal go to the directory you downloaded and run  
`mamba env create --file workflow/envs/snakemake.yaml` (use `conda` instead of `mamba` if you installed Miniconda).

The environments for particular software used by the pipeline will be installed by Snakemake when you run it, so you don't need to install them. The programs in each environment are described in the table below.  
<details>
<summary>Software in the environments used in the pipeline </summary> 

|Environment files | Software | 
| ----: |----: |
|`workflow/envs/snakemake.yaml`|[Snakemake](https://snakemake.github.io/),[Python](https://www.python.org/), [Pandas](https://pandas.pydata.org/), [Click](https://click.palletsprojects.com/en/8.1.x/), [Biopython](https://biopython.org/), [Xonsh](https://xon.sh/), [Graphviz](https://graphviz.org/),[Seqkit](https://bioinf.shenwei.me/seqkit/)|
|`workflow/envs/snippy.yaml`|[Snippy](https://github.com/tseemann/snippy),[Samtools](https://www.htslib.org/)|
|`workflow/envs/liftoff.yaml`|[Litoff](https://github.com/agshumate/Liftoff),[Minimap2]()|
|`workflow/envs/agat.yaml`|[AGAT](https://github.com/NBISweden/AGAT),[Seqkit](https://bioinf.shenwei.me/seqkit/)|
|`workflow/envs/samtools.yaml`|[Samtools](https://www.htslib.org/), [Bedtools](https://bedtools.readthedocs.io/en/latest/index.html), [Bcftools](https://samtools.github.io/bcftools/bcftools.html),[Xonsh](https://xon.sh/),[Pandas](https://pandas.pydata.org/), [Click](https://click.palletsprojects.com/en/8.1.x/), [SciPy](https://scipy.org/), [NumPy](https://numpy.org/) |
|`workflow/envs/depth.yaml`|[Mosdepth](https://github.com/brentp/mosdepth)|
|`workflow/envs/repeatmasker.yaml`|[RepeatMasker](https://www.repeatmasker.org/),[RepeatModeler](https://www.repeatmasker.org/RepeatModeler/), [Bedtools](https://bedtools.readthedocs.io/en/latest/index.html), [Seqkit](https://bioinf.shenwei.me/seqkit/)|
|`workflow/envs/r.yaml` | R, tidyverse, svglite, scales, RColorBrewer||
|`workflow/envs/variants.yaml`| [SnpEff](https://pcingola.github.io/SnpEff/),[DuckDB](https://duckdb.org/), [PyVCF](https://pyvcf.readthedocs.io/en/latest/), [Xonsh](https://xon.sh/),[Pandas](https://pandas.pydata.org/), [Click](https://click.palletsprojects.com/en/8.1.x/), [Biopython](https://biopython.org/), [Bedtools](https://bedtools.readthedocs.io/en/latest/index.html), [Bcftools](https://samtools.github.io/bcftools/bcftools.html)|
</details>

## Configuration

### Input files

#### Tables

* `config/metadata.csv`: A comma-separated table with one sample per row.  
Mandatory columns:  
  * `sample`: sample ID used in the FASTQ file names (no special characters or spaces)  
  * `lineage`: lineage or group name that associates the sample with a reference genome (no special characters or spaces)  
  * `strain`: strain name (a "common name" for each sample, it can be the same as `sample` if you don't have a different one).   
If the plotting will be activated you need one metadata column to color the samples, specify the name of this column in  `metadata2color` in `config/config.yaml`. More columns with free format are allowed. [Example](https://github.com/magwenelab/DiversityPipeline/blob/main/config/metadata.csv). Specify the path to this file in `metadata` in `config/config.yaml`, the default is `config/metadata.csv`.  

* `config/chromosomes.csv`: A comma-separated table with one row per chromosome per lineage.  
Mandatory columns:  
`lineage`: Lineage name (the same as in the metadata table and the names of the reference files)  
`accession`: Sequence ID of the chromosomes in the FASTA and GFF of the reference of each lineage. Make sure each chromosome ID is not repeated in this file.     
`chromosome`: Common name of the chromosome, e.g. chr01, 1, VNI_chr01.  
[Example](https://github.com/magwenelab/DiversityPipeline/blob/main/config/chromosomes.csv).
Specify the path to this file in `chromosomes` in `config/config.yaml`, the default is `config/chromosomes.csv`.

* `config/RepBase.fasta`: Database of repetitive sequences to use for RepeatModeler and RepeatMasker in FASTA format. Needed if the CNV, plotting or database modules with be activated. Specify the path tho this file in `cnv: repeats_database` in `config/config.yaml`.  

* `config/loci.csv`: If you want gene features to be plotted to the depth and MAPQ plots, provide a CSV with the first column `gene_id` with the gene IDs, and the second column `feature` with the name of the feature (locus, pathway, centromere, individual gene name etc.) the gene belongs to. [Example](https://github.com/magwenelab/DiversityPipeline/blob/main/config/loci.csv). Specify the path to this file in `plotting: loci` in `config/config.yaml`.  

* `config/exclude.txt`: If you want to exclude from all analysis some of the samples in your metadata file you can provide a file with a list of sample names to exclude. Specify the path to this file in `samples_to_exclude:` in `config/config.yaml` (no default).

#### Data

* FASTQ files: Paired-end short-read FASTQ files, one forward and one reverse file for each sample. The names of these files should be the names used in the metadata `sample` column, followed by an extension specified in `fastqs: fastq_suffix1` and  `fastqs: fastq_suffix1` in `config/config.yaml`. This files can be gzip compressed. The FASTQ files for all samples should be in the same directory. Specify the path to this directory in `fastqs: directory` (default is `data/samples`).  

* Reference genomes:    
  * If you will use reference genomes with annotation: Provide the FASTA and GFF files for each reference genome. The names of the files must be the ones in the `lineage` column of the metadata (e.g. `VNI.fasta` and `VNI.gff`). Put all the files in the same directory. Specify the path to this directory in `references: directory` in `config/config.yaml` (default is `data/references`).   
  * If you will use a main reference to annotate the reference genomes: Provide the FASTA file for each reference genome. The names of the files must be the ones in the `lineage` column of the metadata, e.g. `VNI.fasta`. Put all the files in the same directory. Specify the path to this directory in `references: directory` in `config/config.yaml` (default is `data/references`). And provide a FASTA and GFF files for the main reference, put both files in the same directory. Specify the path to it in `annotate_references: directory` (default is `data/main_reference`) and specify the names of the files in `annotate_references:fasta` and `annotate_references:gff` (no defaults) in the `config/config.yaml`.   


### Edit the configuration file

Edit the provided `config/config.yaml` file to:
* Select the workflow to run. The `analysis` workflow will run the analysis for one dataset. If you have the complete results (database module activated) of the analysis workflow for two datasets you can join them to create a database with both of them using the `join_datasets` workflow.  
* Provide the paths to the input files and output directory.  
* Activate each module and specify its parameters. See the description of the output below to know which files are created by each module. Activating the `database` module automatically activates the modules `cnv`, `genes_mapq_depth` and `snpeff`.  

## Execution
In a terminal the working directory must be the directory you downloaded.  
Activate the Snakemake environment: `conda activate snakemake`.  
Run the pipeline: `snakemake --cores <n> --sdm conda -p`
  * Snakemake options:  
    * `--cores <n>`: Number of cores you want to use. Mandatory.
    * `--sdm conda`: Specify that Snakemake will use conda environments to run the rules. Mandatory.
    * `-p`: Print the command-lines run by each job in the standard output. Optional.
    * `--conda-frontend conda`: Use it if you are using conda instead of mamba. Mamba is the default.
    * `--rerun-incomplete`: Use it when a past run of the workflow was aborted and you are repeating the run.
    * `--keep-going`: Use it to avoid stoping the workflow when a job fails (i.e. run everything that can run).
    * `-n`: Dry run.  
  


## Output

### Processing of reference genomes

| Path | Description |
| :---------------- | ----: |
| 03.References/{lineage}.gff | Polished GFF file of the reference genome, result of Liftoff annotation using the main reference. Positions are 1-Based. ||
| 03.References/{lineage}_repeats.bed | BED file of regions with repetitive sequences identified by RepeatMasker. Each region is the intersection of diferent types of repetitive sequences identified. Positions are 0-Based. Columns are Accession, Start, End, Types (comma separated list of types in the region).

<details>
<summary> Intermediate files </summary> 

| Path | Description |
| :---------------- | ----: |
| 04.Intermediate_files/03.References/agat_config.yaml | Config file for AGAT
| 04.Intermediate_files/03.References/all_lineages.gff.tsv | Table version of GFF files of all the reference genomes concatenated. Positions are 1-Based. Features are not sorted like in the GFF. Column names different from standard GFF:         accession (seq_id), feature_id (ID), gene_name (Name), gene_id (locus), old_feature_id (original ID before fixing), lineage, and matches_ref_protein,	missing_start_codon,	missing_stop_codon,	inframe_stop_codon (see Liftoff output for last 4).|
| 04.Intermediate_files/03.References/chromosomes.csv | Same as chromosome.csv provided by the user.|
| 04.Intermediate_files/03.References/{main_reference}_fixed_description.gff | GFF with description tag instead of product tag|
| 04.Intermediate_files/03.References/{main_reference}_fixed_ID.gff | GFF with fixed IDs |
| 04.Intermediate_files/03.References/{main_reference}_fixed_locus.gff | GFF with locus tag added |
| 04.Intermediate_files/03.References/{main_reference}_fixed.tsv | Table version of fixed_description GFF |
| 04.Intermediate_files/03.References/{main_reference}.gff | Final fixed GFF with new IDs in the shape of `<locus>-<level2 tag_level2 number>-<level3 tag and number>` |
| 04.Intermediate_files/03.References/{main_reference}.tsv | TSV version of fixed GFF |
| 04.Intermediate_files/03.References/{lineage}/{main_reference}.fasta | Symlink to original FASTA |
| 04.Intermediate_files/03.References/{lineage}/{main_reference}.fasta.fai | FASTA index created by Liftoff |
| 04.Intermediate_files/03.References/{lineage}/{main_reference}.gff | Symlink to fixed GFF |
| 04.Intermediate_files/03.References/{lineage}/{main_reference}.gff_db | DB of GFF created by Liftoff |
| 04.Intermediate_files/03.References/{lineage}/intermediate_liftoff/ | See [Liftoff output](https://github.com/agshumate/Liftoff?tab=readme-ov-file#usage) |
| 04.Intermediate_files/03.References/{lineage}/liftoff.gff | GFF from Liftoff before polishing |
| 04.Intermediate_files/03.References/{lineage}/unmapped_features.txt | List of features not found in reference genome |
| 04.Intermediate_files/03.References/{lineage}/{lineage}.cds.fa | Nucleotide sequences of all isoforms in reference genome |
| 04.Intermediate_files/03.References/{lineage}/{lineage}.fasta.fai 
| 04.Intermediate_files/03.References/{lineage}/{lineage}.fasta.index
| 04.Intermediate_files/03.References/{lineage}/{lineage}.fasta.mmi
| 04.Intermediate_files/03.References/{lineage}/{lineage}.gff.tsv | Table version of polished GFF of reference genome. Positions are 1-Based. Features are not sorted like in the GFF. |
| 04.Intermediate_files/03.References/{lineage}/{lineage}.prots.fa | Protein sequences of all isoforms in reference genome |
| 04.Intermediate_files/03.References/{lineage}/repeats/01_simple/{lineage}.bed | BED file of simple repetitive sequences. Positions are 0-Based. |
| 04.Intermediate_files/03.References/{lineage}/repeats/01_simple/ | See [RepeatMasker output](https://www.repeatmasker.org/webrepeatmaskerhelp.html)|
| 04.Intermediate_files/03.References/{lineage}/repeats/02_complex/{lineage}.bed | BED file of complex repetitive sequences. Positions are 0-Based. |
| 04.Intermediate_files/03.References/{lineage}/repeats/02_complex/ | See [RepeatMasker output](https://www.repeatmasker.org/webrepeatmaskerhelp.html)|
| 04.Intermediate_files/03.References/{lineage}/repeats/03_known/{lineage}.bed| BED file of known repetitive sequences. Positions are 0-Based. |
| 04.Intermediate_files/03.References/{lineage}/repeats/03_known/ | See [RepeatMasker output](https://www.repeatmasker.org/webrepeatmaskerhelp.html)|
| 04.Intermediate_files/03.References/{lineage}/repeats/04_unknown/{lineage}.bed | BED file of unknown repetitive sequences. Positions are 0-Based. |
| 04.Intermediate_files/03.References/{lineage}/repeats/04_unknown/ | See [RepeatMasker output](https://www.repeatmasker.org/webrepeatmaskerhelp.html)|
| 04.Intermediate_files/03.References/{lineage}/repeats/{lineage}_db/ | Database created with RepeatModeler's BuildDatabase |
| 04.Intermediate_files/03.References/{lineage}/repeats/{lineage}_known.fa | FASTA file of known families of repetitive sequences identified by RepeatModeler |
| 04.Intermediate_files/03.References/{lineage}/repeats/{lineage}_unknown.fa | FASTA file of unknown families of repetitive sequences identified by RepeatModeler |
| 04.Intermediate_files/03.References/{lineage}/repeats/RModeler/ | See [RepeatModeler output](https://www.repeatmasker.org/RepeatModeler/)|

</details>

### Snippy

| Path | Description |
| :---------------- | ----: |
| 01.Samples/snippy/{sample}/snps.bam | BAM file of alignment between short reads of sample and corresponding reference genome. |
| 01.Samples/snippy/{sample}/snps.consensus.fa | FASTA file of the reference genome with all variants instantiated. |
| 01.Samples/snippy/{sample}/snps.vcf | Called variants in VCF format. Positions are 01-Based.|
| 01.Samples/snippy/{sample}/* | Other files from the [Snippy output](https://github.com/tseemann/snippy?tab=readme-ov-file#output-files).|


### Depth and quality

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 01.Samples/depth_quality/{sample}/depth_by_chrom_good.tsv
| 01.Samples/depth_quality/{sample}/depth_by_chrom_raw.tsv
| 01.Samples/depth_quality/{sample}/mapping_stats.tsv
| 01.Samples/depth_quality/{sample}/mapq_depth_window.bed | Positions are 0-Based.
| 02.Dataset/depth_quality/depth_by_chrom_good.tsv
| 02.Dataset/depth_quality/depth_by_chrom_raw.tsv
| 02.Dataset/depth_quality/mapping_stats.tsv

<details>
<summary> Intermediate files </summary> 

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage_good.mosdepth.global.dist.txt
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage_good.mosdepth.region.dist.txt
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage_good.mosdepth.summary.txt
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage_good.regions.bed.gz
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage_good.regions.bed.gz.csi
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage.mosdepth.global.dist.txt
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage.mosdepth.region.dist.txt
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage.mosdepth.summary.txt
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage.regions.bed.gz
| 04.Intermediate_files/01.Samples/mosdepth/{sample}/coverage.regions.bed.gz.csi
| 04.Intermediate_files/01.Samples/depth_quality/{sample}/depth_by_windows.tsv
| 04.Intermediate_files/01.Samples/depth_quality/{sample}/depth_distribution.tsv
| 04.Intermediate_files/01.Samples/depth_quality/{sample}/mapq.bed
| 04.Intermediate_files/01.Samples/depth_quality/{sample}/mapq_window.bed
| 04.Intermediate_files/01.Samples/depth_quality/{sample}/snps_good.bam
| 04.Intermediate_files/01.Samples/depth_quality/{sample}/snps_good.bam.bai
| 04.Intermediate_files/01.Samples/filtered_samples/{sample}.txt
| 04.Intermediate_files/02.Dataset/depth_quality/unfiltered_mapping_stats.tsv
| 04.Intermediate_files/03.References/filtered_lineages/{lineage}.txt

</details>

### Annotation

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 01.Samples/annotation/{sample}/annotation.gff | Positions are 1-Based.
| 01.Samples/annotation/{sample}/cds.fa
| 01.Samples/annotation/{sample}/proteins.fa

<details>
<summary> Intermediate files </summary> 

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 04.Intermediate_files/01.Samples/liftoff/{sample}/intermediate_liftoff/(see liftoff output)
| 04.Intermediate_files/01.Samples/liftoff/{sample}/lifted.gff
| 04.Intermediate_files/01.Samples/liftoff/{sample}/ref.gff_db
| 04.Intermediate_files/01.Samples/liftoff/{sample}/unmapped_features.txt
| 04.Intermediate_files/01.Samples/annotation/{sample}/cds.csv
| 04.Intermediate_files/01.Samples/annotation/{sample}/proteins.csv
| 04.Intermediate_files/02.Dataset/sequences.csv
| 04.Intermediate_files/agat_config.yaml

</details>

### Depth and quality of genes

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 01.Samples/depth_quality/{sample}/feature_mapq_depth.tsv
| 02.Dataset/depth_quality/feature_mapq_depth.tsv

### CNV calling

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 01.Samples/cnv/{sample}/cnv_calls.tsv | Positions are 0-Based.|
| 02.Dataset/cnv/cnv_calls.tsv | Positions are 0-Based.|

### SNP effects

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 02.Dataset/snps/effects.tsv
| 02.Dataset/snps/lofs.tsv
| 02.Dataset/snps/nmds.tsv
| 02.Dataset/snps/presence.tsv
| 02.Dataset/snps/variants.tsv | Positions are 1-Based.

<details>
<summary> Intermediate files </summary> 

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 04.Intermediate_files/02.Dataset/snps/{lineage}_effects.tsv
| 04.Intermediate_files/02.Dataset/snps/{lineage}_intersection.vcf | Positions are 1-Based.
| 04.Intermediate_files/02.Dataset/snps/{lineage}_lofs.tsv
| 04.Intermediate_files/02.Dataset/snps/{lineage}_nmds.tsv
| 04.Intermediate_files/02.Dataset/snps/{lineage}_presence.tsv
| 04.Intermediate_files/02.Dataset/snps/{lineage}_snpeff.genes.txt
| 04.Intermediate_files/02.Dataset/snps/{lineage}_snpeff.html
| 04.Intermediate_files/02.Dataset/snps/{lineage}_snpeff.vcf
| 04.Intermediate_files/02.Dataset/snps/{lineage}_variants.tsv | Positions are 1-Based.
| 04.Intermediate_files/03.References/snpeff_data/Cryptococcus_neoformans_{lineage}/
| 04.Intermediate_files/03.References/snpeff_data/{lineage}.done
| 04.Intermediate_files/03.References/snpeff_data/snpEff.config

</details>

### Database

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 02.Dataset/database.db

### Plots

| Path | Description |
| :---------------- | ----: |
| 01.Samples/plots/{sample}/depth_by_chrom.png
| 01.Samples/plots/{sample}/depth_by_windows.png
| 01.Samples/plots/{sample}/depth_chrom_distribution.png
| 01.Samples/plots/{sample}/depth_global_distribution.png
| 01.Samples/plots/{sample}/mapq.png
| 02.Dataset/plots/dataset_depth_by_chrom.png
| 02.Dataset/plots/dataset_summary.png

<details>
<summary> Intermediate files </summary> 

| Path | Description | Column names |
| :---------------- | ----: |----: |
| 04.Intermediate_files/03.References/loci.csv
| 04.Intermediate_files/03.References/loci_to_plot.tsv | Positions are 1-Based.

</details>

## Filegraph


![Filegraph](all.svg)
