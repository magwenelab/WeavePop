
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

  * `config.yaml`: This is a Snakemake configuration file, provided [here](https://github.com/magwenelab/DiversityPipeline/blob/main/config.yaml), that you must edit according to your filenames and desired parameters.
  * `config/sample_metadata.csv` with columns (using these exact names):  sample (the names in the fastq file names), group (lineage or group to associate to a reference genome), strain, more-optional-metadata-fields.
  * Optional: `config/loci.csv`, CSV with first column with the gene IDs and second column with the name of the locus/pathway the gene belongs to if you want genes to be plotted to coverage and MAPQ plots  No column names.
  * `config/chromosome_names.csv` with columns (without column names): group, chromosome ID (the sequence ID in the Fasta and GFF of the references), chromosome name (typically a number). If your genomes are Complete Genomes from NCBI use `bash get-chromosome_names.sh` to get this file.
  * `config/features.txt` list of feature names to lift over. This file is provided in this repository.
  * `references/` directory with:
    * Fasta files to use as reference genomes. The names of the files must be the ones in the "group" column of the `files/sample_metadata.csv`, e.g. `VNI.fasta`
    * Optional: Fasta and GFF files of main reference (one with available annotation with the desired gene IDs). If your main reference is one of your reference genomes, duplicate the genome files and give another name to the ones that will be used as main.
    * If main reference is not provided, GFF files of reference genomes are needed (with the same name as the fastas).
      * If your genomes have a mitochondrial chromosome you can run `bash get-removed-chrom.sh path-to-fasta path-to-gff seq_id` to remove it, in an environment with Seqkit available (diversity).

See example files [here](https://github.com/magwenelab/DiversityPipeline/tree/main/files).

