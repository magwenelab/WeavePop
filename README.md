
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

To install the required programs for each module you can run `mamba env create --file envs/envname.yaml` using the specified environment file in the table bellow. The environments for Modules 1 - 3 are installed by Snakemake, so you don't need to do it.

| Module | Software | Environment file |
| :---------------- | ----: |----: |
| Module 0|[Sra-Tools](https://github.com/ncbi/sra-tools) , [Entrez-Direct](https://www.ncbi.nlm.nih.gov/books/NBK25501/) |`envs/sra-tools.yaml`|
| Module 1|[Litoff](https://github.com/agshumate/Liftoff), [AGAT](https://github.com/NBISweden/AGAT)|`envs/liftoff.yaml`,`envs/agat.yaml`|
| Module 2|[Snippy](https://github.com/tseemann/snippy), [Litoff](https://github.com/agshumate/Liftoff), [AGAT](https://github.com/NBISweden/AGAT)|`envs/snippy.yaml`, `envs/liftoff.yaml`, `envs/agat.yaml`|
| Module 3|[Mosdepth](https://github.com/brentp/mosdepth), [Samtools](https://www.htslib.org/)|`envs/depth.yaml`|
| Module 3|R and R libraries -- tidyverse ComplexHeatmap, svglite, scales, RColorBrewer|`envs/r.yaml`|


## Structure of repository:  
  * The working directory has the scripts and Snakefiles to run.  
  * `files/` has some of the starting files and files created by the pipeline.
  * `scripts/` has the scripts used by the Snakefiles, not by the user directly.  
  * `references/` has the reference genomes.  
  * `analyses/` has one directory per sample, all the resulting files of the analyses performed per sample will be there with a generic name.  
  * `results/` has the resulting files of the analyses that consider all the samples.  
  * `logs/` has the log files of all runs.  

## Starting files: 
  * `files/sample_metadata.csv` with columns (using these names):  sample (the names in the fastq file names), group (lineage or group to associate to a reference genome), strain, more-optional-metadata-fields.
  * Lists of genes of loci of interest:  
    * `files/locusA.txt` (with IDs of genes in main reference GFF)
  * `files/chromosome_names.csv` with columns (without column names): group, chromosome ID (the sequence ID in the Fasta and GFF of the references), chromosome name (typically a number). If your genomes are Complete Genomes from NCBI use `bash get-chromosome_names.sh` to get this file.
  * `files/features.txt` list of feature names to lift over. This file is provided in this repository.
  * `references/` directory with:
    * Fasta files to use as reference genomes. The names of the files must be the ones in the "group" column of the `files/sample_metadata.csv`, e.g. `VNI.fasta`
    * Optional: Fasta and GFF files of main reference (one with available annotation with the desired gene IDs). If your main reference is one of your reference genomes, duplicate the genome files and give another name to the ones that will be used as main.
    * If main reference is not provided, GFF files of reference genomes are needed (with the same name as the fastas).
      * If your genomes have a mitochondrial chromosome you can run `bash get-removed-chrom.sh path-to-fasta path-to-gff seq_id` to remove it, in an environment with Seqkit available (crypto_div).


## Modules to be run in this order:
You can start from Module 0, Module 1 or Module 2.

### Module 0 (Optional):
 To download all FASTQs of a BioProject, compress them and rename them with the SRS sample ID. If a sample has more than one set of read files they will be concatenated into one set (one forward and one reverse fastq).

| Input provenance | Input | Description |
| :---------------- | ----: | ----:|
| - | BioProject ID |BioProject with valid SRA short read sequencing fastq files.|

| Output | Output description |Needed for:|
| :---------------- | ----: |----: | 
|`srafiles/`|Directory with `.sra` files.| This module.|
|`fastqs/`|Directory with `.fastq` files|This module.|
|`files/read_pair_table.csv`|Column names: sample, run, file1, file2, size|This module.|
|`files/unpaired_fastqs.csv`|Column names: sample, run, files|This module.|
|`files/samples.txt`|List of SRS sample IDs|This module.|
| `fastq_combined/` |Directory with `.fq.gz` files|Modules 1-3|


~~~
$ conda activate sra-tools
$ xonsh get-seqdata-of-bioproject.xsh -p PRJNA685103
$ parallel xonsh get-fastqs-combined.xsh {} files/read_pair_table.csv fastqs/ fastq_combined/ :::: files/samples.txt 
~~~

It is possible that this module does not download all the samples that you expect (sometimes a BioProject has a BioSample associated but there is no SRA file for it). Make sure the **metadata** you provide only has the samples that will be processed in Module 2.

### Module 1 (Optional):
Lift over annotations from the main reference into the reference genomes.  

| Input provenance | Input | Description |
| :---------------- | ----: | ----:|
||||

| Output | Description |Needed for:|
| :---------------- | ----: |----: | 
| `references/{lineage}.gff`|||
|`references/{lineage}.gff.tsv`|||
|`references/references_unmapped.svg`|||

~~~
$ conda activate diveristy
$ snakemake --snakefile Snakefile-references.smk --cores 1 --use-conda -p 
~~~

⚠️ `--cores 1` is because there is a problem if Liftoff runs in parallel because the different jobs try to create `mainReference.gff_db` at the same time and that is not cool.     

### Module 2:
Run the analysis per sample. It runs **snippy**, **liftoff** and **agat** for each sample, it **extracts sequences** (cds and protein) of each sample and **concatenates** them by cds and by protein.

| Input | Description |Input provenance |
|:---- | ----: |----------------:|
||||

| Output | Description |Needed for:|
| :---------------- | ----: |----: | 
|`{lineage}.gff.tsv`|||
|`{lineage}_predicted_cds.fa`|||
|`{lineage}_predicted_proteins.fa`|||
|`analyses/{sample}/snps.consensus.fa` and extra assembly files|||
|`analyses/{sample}/snps.bam` and extra alignment files |||
|`analyses/{sample}/snps.vcf` and extra variant calling files |||
|`analyses/{sample}/lifted.gff_polished` and extra annotation files|||
|`analyses/{sample}/predicted_cds.fa`  and extra index files|||
|`analyses/{sample}/predicted_proteins.fa`  and extra index files|||
|`results/cds/{protein}.fa`|||
|`results/proteins/{protein}.fa`|||
~~~
$ conda activate diveristy
$ snakemake --snakefile Snakefile-main.smk --cores <n> --use-conda -p 
~~~

### Module 3: Quality and depth analyses

Generates **quality and coverage** plots, and adds the MAPQ and Coverage of each locus's window to the GFF file.  
| Input | Description |Input provenance |
|:---- | ----: |----------------:|
||||

| Output | Description |Needed for: |
|:---- | ----: |----------------:|
|`results/mapped_reads.svg` and `results/mapping_stats.txt`|plot and table with fraction of mappied reads per sample.||
|`analyses/{sample}/snps.bam.stats` and `analyses/{sample}/bamstats/`|directory with `plot-bamstats` resulting plots.  ||
|`analyses/{sample}/coverage.svg`|coverage along chromosome plot (with location of interesting loci).  ||
|`analyses/{sample}/cov_distribution_.svg`| ditribution of coverage values plot.  ||
|`analyses/{sample}/mapq.svg` |mapping quality along chromosome plot (with location of interesting loci).  ||
|`analyses/{sample}/mapq_distribution.svg`| distribution of maping quality values plot.||
|`analyses/{sample}/annotation.gff` |GFF file with complete annotation plus average MAPQ and coverage of windows in which the features are located. |||
|`results/cov_norm_good.csv` |table with all coverage stats of good quality mappings, per chromosome of all samples (genome-wide and per chromosome mean and median and normalized) ||
|`results/cov_global_good.svg` |plot of mean and median genome-wide coverage of good quality mappings of all samples.||
|`results/cov_median_good.svg`| plot of median coverage per chromosome (of good quality mappings) normalized by genome-wide median coverage.||
|`results/cov_mean_good.svg` |plot of mean coverage per chromosome (of good quality mappings) normalized by genome-wide mean coverage.||
|`results/cov_norm_raw.csv`| table with all coverage stats of all mappings, per chromosome of all samples (genome-wide and per chromosome mean and median and normalized) ||
|`results/cov_global_raw.svg` |plot of mean and median genome-wide coverage of all mappings of all samples.||
|`results/cov_median_raw.svg`| plot of median coverage per chromosome (of all mappings) normalized by genome-wide median coverage.||
|`results/cov_mean_raw.svg`| plot of mean coverage per chromosome (of all mappings) normalized by genome-wide mean coverage.||
 
~~~
$ conda activate diveristy
$ snakemake --snakefile Snakefile-depth-quality.smk --cores <n> --use-conda -p 
~~~

# To run all

```
echo "Running References" &&
snakemake --snakefile Snakefile-references.smk --cores 1 --use-conda --conda-frontend conda -p &> all.log &&
echo "Running Main" &&
snakemake --snakefile Snakefile-main.smk --cores 12 --use-conda --conda-frontend conda -p  &>> all.log &&
echo "Running Depth" &&
snakemake --snakefile Snakefile-depth-quality.smk --cores 8 --use-conda --conda-frontend conda -p &>> all.log
```
