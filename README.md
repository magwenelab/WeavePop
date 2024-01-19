
# Reference-based mapping, variant calling and assembly / Annotation liftover / Mapping quality and coverage

## Broad description


## Requirements

The environment from which Modules 1-3 must me run has the following software and you can install it with: `mamba env create --file envs/diversity.yml`
* Mamba/Conda [Microforge3](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
* Python
* Python modules -- Pandas, Click
* [Xonsh](https://xon.sh/)
* [Snakemake](https://snakemake.github.io/)
* [Graphviz](https://graphviz.org/) (optional, to seethe Snakemake DAG in a graph) 
* [Seqkit](https://bioinf.shenwei.me/seqkit/)

The environments for Modules 1 - 3 are installed by Snakemake, so you don't need to do it. If you want to install them you can run `mamba env create --file envs/envname.yaml` using the specified environment file in the table bellow. 
<details>
<summary>Requirements per module </summary> 

| Module | Software | Environment files |
| :---------------- | ----: |----: |
| Module 0|[Sra-Tools](https://github.com/ncbi/sra-tools) , [Entrez-Direct](https://www.ncbi.nlm.nih.gov/books/NBK25501/) |`envs/sra-tools.yaml`|
| Module 1|[Litoff](https://github.com/agshumate/Liftoff), [AGAT](https://github.com/NBISweden/AGAT)|`envs/liftoff.yaml`,`envs/agat.yaml`|
| Module 2|[Snippy](https://github.com/tseemann/snippy), [Litoff](https://github.com/agshumate/Liftoff), [AGAT](https://github.com/NBISweden/AGAT)|`envs/snippy.yaml`, `envs/liftoff.yaml`, `envs/agat.yaml`|
| Module 3|[Mosdepth](https://github.com/brentp/mosdepth), [Samtools](https://www.htslib.org/)|`envs/depth.yaml`, `envs/samtools.yaml`|
| Module 4|[Samtools](https://www.htslib.org/), Gnuplot, matplotlib, tectonic, texlive-core, R and R libraries -- tidyverse ComplexHeatmap, svglite, scales, RColorBrewer|`envs/plot-bamstats.yaml`,`envs/r.yaml`|
</details>

## Structure of repository:  
  * The working directory has the scripts and Snakefiles to run.  
  * `files/` has some of the starting files and files created by the pipeline.
  * `scripts/` has the scripts used by the Snakefiles, not by the user directly.  
  * `references/` has the reference genomes.  
  * `analyses/` has one directory per sample, all the resulting files of the analyses performed per sample will be there with a generic name.  
  * `results/` has the resulting files of the analyses that consider all the samples.  
  * `logs/` has the log files of all runs.  

## Starting files: 
  * `config.yaml`: This is a Snakemake configuration file, provided in this repository, that you must edit according to your filenames and desired parameters.
  * `files/sample_metadata.csv` with columns (using these names):  sample (the names in the fastq file names), group (lineage or group to associate to a reference genome), strain, more-optional-metadata-fields.
  * Lists of genes of loci of interest:  
    * `files/locusA.txt` (with IDs of genes in main reference GFF)
  * `files/chromosome_names.csv` with columns (without column names): group, chromosome ID (the sequence ID in the Fasta and GFF of the references), chromosome name (typically a number). If your genomes are Complete Genomes from NCBI use `bash get-chromosome_names.sh` to get this file.
  * `files/features.txt` list of feature names to lift over. This file is provided in this repository.
  * `references/` directory with:
    * Fasta files to use as reference genomes. The names of the files must be the ones in the "group" column of the `files/sample_metadata.csv`, e.g. `VNI.fasta`
    * Optional: Fasta and GFF files of main reference (one with available annotation with the desired gene IDs). If your main reference is one of your reference genomes, duplicate the genome files and give another name to the ones that will be used as main.
    * If main reference is not provided, GFF files of reference genomes are needed (with the same name as the fastas).
      * If your genomes have a mitochondrial chromosome you can run `bash get-removed-chrom.sh path-to-fasta path-to-gff seq_id` to remove it, in an environment with Seqkit available (diversity).


## Modules
Modules should be run in order 0 to 4.  
You can start from Module 0, Module 1 or Module 2.

### Module 0 (Optional): Download all FASTQs of a BioProject
This module has one script (`get-seqdata-of-bioproject.xsh`) that takes a BioProject ID uses Entrez-Direct to know which samples from this project are in the [SRA](https://www.ncbi.nlm.nih.gov/sra/docs/), and Sra-Tools to download them in FASTQ format. 
The second script (`get-fastqs-combined.xsh`) renames and compresses the FASTQ files, if a sample has more than one set of read files they will be concatenated into one set (one forward and one reverse FASTQ). If there are single-end sequencing reads they will be ignored.

<details>
<summary>Input</summary> 

| Input | Description |
| ----: | ----:|
| BioProject ID |BioProject with valid short read paired-end sequencing FASTQ files in the SRA|
</details>
<details>
<summary>Output</summary> 

| Output | Output description |Needed for|
| :---------------- | ----: |----: | 
|`srafiles/`|Directory with `.sra` files.| This module|
|`fastqs/`|Directory with `.fastq` files|This module|
|`files/read_pair_table.csv`|Column names: sample, run, file1, file2, size|This module|
|`files/unpaired_fastqs.csv`|Column names: sample, run, files|This module|
|`files/samples.txt`|List of SRS sample IDs|This module|
|`fastq_combined/` |Directory with `.fq.gz` files|Modules 2|
</details>

Run:
~~~
$ conda activate sra-tools
$ xonsh get-seqdata-of-bioproject.xsh -p PRJNA685103
$ parallel xonsh get-fastqs-combined.xsh {} files/read_pair_table.csv fastqs/ fastq_combined/ :::: files/samples.txt 
~~~

It is possible that this module does not download all the samples that you expect, sometimes a BioProject has a BioSample associated but there is no SRA file for it.

### Module 1 (Optional): Annotate reference genomes
Lift over annotations from a main reference (a genome with available FASTA and GFF files and the gene IDs that you want to use) into the reference genomes, of one or more group/lineage, that will be used for mapping the reads.  
<details>
<summary>Input</summary> 

| Input | Description |Input origin |
| :---- | ----:|----------------: |
|`files/sample_metadata.csv`| Columns: `sample`(the names in the fastq file names), `group` (lineage or group to associate to a reference genome), `strain`, more-optional-metadata-fields|You|
|`references/mainReference.fasta`|Main reference genome assembly.|You (Tipically public genome from FungiDB or NCBI)|
|`references/mainReference.gff`|GFF annotation file of the genome described above|You (Tipically public genome from FungiDB or NCBI)|
|`references/{lineage}.fasta`|Genome assembly of each group/lineage. The filename must have the names used in the `group` column of the `sample_metadata.csv`|You|
|`files/features.txt`|List of level 1 features to lift over from the references, check [this](https://github.com/agshumate/Liftoff?tab=readme-ov-file#feature-types) to know more |Provided in this repository|
</details>

<details>
<summary>Output</summary> 

|Output | Description |Needed for|
|:---------------- | ----: |----: | 
|`references/mainReference.gff.tsv`|Tabular version of the GFF of the main reference|Module 2,3|
|`references/{lineage}.gff`|Annotation GFF of each group/lineage|Module 2,3|
|`references/references_unmapped.svg`|Heatmap showing the features that were not lifted over from the main reference to each group's genome|-|
</details>

Run:
~~~
$ conda activate diversity
$ snakemake --snakefile Snakefile-references.smk --cores 1 --use-conda -p 
~~~

⚠️ `--cores 1` is because there is a problem if Liftoff runs in parallel because the different jobs try to create `mainReference.gff_db` at the same time and that is not cool.     

### Module 2: Mapping, annotation and sequence extraction
This module runs Snippy to map the reads of each sample to the reference genome of the corresponding group/lineage, call variants and provide a genome assembly. It uses Liftover to do a functional annotation, and it extracts the nucleotide and aminoacid sequences of all the genes with AGAT.

<details>
<summary>Input</summary>  

| Input | Description |Input origin |
|:---- | ----: |----------------:|
|`files/sample_metadata.csv`| Columns: `sample`(the names in the fastq file names), `group` (lineage or group to associate to a reference genome), `strain`, more-optional-metadata-fields|You|
|`fastq_combined/` |Directory with `.fq.gz` files|Module 0 or You |
|`references/{lineage}.fasta`|Genome assembly of each group/lineage. The filename must have the names used in the `group` column of the `sample_metadata.csv`|You|
|`references/{lineage}.gff`|Annotation GFF of each group/lineage|Module 1 or You|
|`files/features.txt`|List of level 1 features to lift over from the references, check [this](https://github.com/agshumate/Liftoff?tab=readme-ov-file#feature-types) to know more |Provided in this repository|
</details>

<details>
<summary>Output</summary> 

| Output | Description |Needed for|
| :---------------- | ----: |----: | 
|`{lineage}_predicted_cds.fa`|Fasta file with the coding sequences (for each isoform) for each reference genome||
|`{lineage}_predicted_proteins.fa`|Fasta file with aminoacid sequence of each protein isoform for each reference genome||
|`analyses/{sample}/snps.consensus.fa` and extra assembly files|A version of the corresponding reference genome with the SNPs of the sample instead.||
|`analyses/{sample}/snps.bam` and extra alignment files |Mapping BAM file for each sample|Module 3|
|`analyses/{sample}/snps.vcf` and extra variant calling files |Variants VCF file for each sample||
|`analyses/{sample}/lifted.gff_polished` and extra annotation files|GFF annotation file|Module 3|
|`analyses/{sample}/predicted_cds.fa`  and extra index files|Fasta file with the coding sequences (for each isoform) for each sample||
|`analyses/{sample}/predicted_proteins.fa`  and extra index files|Fasta file with aminoacid sequence of each protein isoform for each sample||
|`results/cds/{protein}.fa`|A fasta file for each isoform with the coding sequence of all samples||
|`results/proteins/{protein}.fa`|A fasta file for each isoform with the aminoacid sequence of all samples||
</details>

~~~
$ conda activate diversity
$ snakemake --snakefile Snakefile-main.smk --cores <n> --use-conda -p 
~~~

### Module 3: Quality and coverage analyses

This module generates mapping quality and coverage plots, and adds the MAPQ and Coverage of each locus's window to the GFF file.  

<details>
<summary>Input</summary>  

| Input | Description |Input origin |
|:---- | ----: |----------------:|
|`files/sample_metadata.csv`| Columns: `sample`(the names in the fastq file names), `group` (lineage or group to associate to a reference genome), `strain`, more-optional-metadata-fields|You|
|`analyses/{sample}/snps.bam`||Module 2|
|`analyses/{sample}/lifted.gff_polished`||Module 2|
|`references/{lineage}.gff`||Module 1 or You|
|`files/{Locus}.txt`| One file per locus of interest that you want to show in the plots of coverage and MAPQ along the chromosomes. The file is a list of gene IDs and the name of the file is the one you want to appear in the plots.|You |
</details>
<details>
<summary>Output</summary> 

| Output | Description |Needed for |
|:---- | ----: |----------------:|
|`analyses/{sample}/coverage.regions.bed.gz`|Mean coverage of each window | This module and Module 4|
|`analyses/{sample}/coverage_good.regions.bed.gz`|Mean coverage of each window without bad quality mappings |Module 4|
|`analyses/{sample}/cov.csv`| Number of positions with each coverage value |Module 4|
|`analyses/{sample}/mapq.csv`| Number of positions with each MAPQ value |Module 4|
|`analyses/{sample}/snps.bam.stats`|Results of samtools stats |This module and Module 4|
|`results/mapping_stats.txt`|Stats on number o mapped reads of all samples|Module 4|
|`analysis/{sample}/mapq.bed`|Mapping quality of each position|Module 4|
|`analysis/{sample}/mapq_window.bed`|Mean mapping quality of each window||
|`analysis/{sample}/mapq_cov_window.bed`|Mean MAPQ and Coverage of each window ||
|`analyses/{sample}/annotation.gff` |GFF file with complete annotation plus average MAPQ and coverage of the windows in which the features are located. ||
|`references/{lineage}.gff.tsv`|GFF in table format |This module|
|`files/loci_to_plot.tsv`|GFF in TSV format of all reference genomes, only with genes in loci of interest. Includes a column with the name of the locus |Module 4|
</details>

~~~
$ conda activate diversity
$ snakemake --snakefile Snakefile-depth-quality.smk --cores <n> --use-conda -p 
~~~

### Module 4: Quality and coverage plotting

The module takes the output files of Module 3 to make mapping quality and coverage plots.

<details>
<summary>Input</summary>  

| Input | Description |Input origin |
|:---- | ----: |----------------:|
|`files/sample_metadata.csv`| Columns: `sample`(the names in the fastq file names), `group` (lineage or group to associate to a reference genome), `strain`, more-optional-metadata-fields|You|
|`files/chromosome_names.csv`||You, if your genomes are Complete Genomes from NCBI run `bash get-chromosome_names.sh` to get this file.|
|`analyses/{sample}/coverage.regions.bed.gz`|Mean coverage of each window |Module 3|
|`analyses/{sample}/coverage_good.regions.bed.gz`|Mean coverage of each window without bad quality mappings |Module 3|
|`analyses/{sample}/cov.csv`|Number of positions with each coverage value |Module 3|
|`analyses/{sample}/mapq.csv`|Number of positions with each MAPQ value |Module 3|
|`analyses/{sample}/snps.bam.stats`|Results of samtools stats|Module 3|
|`results/mapping_stats.txt`|Stats on number o mapped reads of all samples|Module 3|
|`analysis/{sample}/mapq_window.bed`|Mean MAPQ and Coverage of each window |Module 3|
|`files/loci_to_plot.tsv`|GFF in TSV format of all reference genomes, only with genes in loci of interest. Includes a column with the name of the locus|Module 3|
</details>

<details>
<summary>Output</summary>  

| Output | Description |
|:---- | ----: |
|`analyses/{sample}/bamstats/`|Directory with `plot-bamstats` resulting plots. |
|`results/mapped_reads.svg` and `results/mapping_stats.txt`|plot and table with fraction of mappied reads per sample.|
|`analyses/{sample}/coverage.svg`|coverage along chromosome plot (with location of interesting loci).  |
|`analyses/{sample}/cov_distribution.svg`| ditribution of coverage values plot.  |
|`analyses/{sample}/mapq.svg` |mapping quality along chromosome plot (with location of interesting loci).  |
|`analyses/{sample}/mapq_distribution.svg`| distribution of maping quality values plot.|
|`results/cov_norm_good.csv` |table with all coverage stats of good quality mappings, per chromosome of all samples (genome-wide and per chromosome mean and median and normalized) |
|`results/cov_global_good.svg` |plot of mean and median genome-wide coverage of good quality mappings of all samples.|
|`results/cov_median_good.svg`| plot of median coverage per chromosome (of good quality mappings) normalized by genome-wide median coverage.|
|`results/cov_mean_good.svg` |plot of mean coverage per chromosome (of good quality mappings) normalized by genome-wide mean coverage.|
|`results/cov_norm_raw.csv`| table with all coverage stats of all mappings, per chromosome of all samples (genome-wide and per chromosome mean and median and normalized) |
|`results/cov_global_raw.svg` |plot of mean and median genome-wide coverage of all mappings of all samples.|
|`results/cov_median_raw.svg`| plot of median coverage per chromosome (of all mappings) normalized by genome-wide median coverage.|
|`results/cov_mean_raw.svg`| plot of mean coverage per chromosome (of all mappings) normalized by genome-wide mean coverage.|
</details>

~~~
$ conda activate diversity
$ snakemake --snakefile Snakefile-plotting.smk --cores <n> --use-conda -p 
~~~