
# Reference-based mapping, variant calling and assembly / Annotation liftover / Mapping quality and coverage

## Description

This is a Snakemake workflow to map short-reads of samples to the reference genome of the corresponding lineage/group. You will get the mapping file, variant calling file, a reference-based assembly, and an annotation GFF. Additionaly you can get analyses of the coverage and the mapping quality, plots with this analyses' results, including detection of structural variants, and a table with the intersection of the detected genetic variants between the samples of the same group.  If you want to have a common naming scheme of your genes and/or you don't have GFF files of your reference genomes you can provide a main reference to lift over the annotation from this one to your reference genomes.

## Requirements

The environment from which the workflow must be run has the following software and you can install it with: `mamba env create --file envs/diversity.yml`
* Mamba/Conda [Microforge3](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
* Python
* Python modules -- Pandas, Click, Biopython
* [Xonsh](https://xon.sh/)
* [Snakemake](https://snakemake.github.io/)
* [Graphviz](https://graphviz.org/) (optional, to seethe Snakemake DAG in a graph) 
* [Seqkit](https://bioinf.shenwei.me/seqkit/)

The environments for particular software used in Snakemake rules are installed by Snakemake, so you don't need to do it. If you want to install them you can run `mamba env create --file workflow/envs/envname.yaml` using the specified environment file in the table bellow. 
<details>
<summary>Requirements per module </summary> 

| Module | Software | Environment files |
| :---------------- | ----: |----: |
| Module: Annotate references|[Litoff](https://github.com/agshumate/Liftoff)|`workflow/envs/liftoff.yaml`|
| Module: Main |[Snippy](https://github.com/tseemann/snippy), [Litoff](https://github.com/agshumate/Liftoff), [AGAT](https://github.com/NBISweden/AGAT)|`workflow/envs/snippy.yaml`, `workflow/envs/liftoff.yaml`, `workflow/envs/agat.yaml`|
| Module: Coverage - Quality|[Mosdepth](https://github.com/brentp/mosdepth), [Samtools](https://www.htslib.org/), R and R libraries -- tidyverse ComplexHeatmap, svglite, scales, RColorBrewer|`workflow/envs/depth.yaml`, `workflow/envs/samtools.yaml`, `workflow/envs/r.yaml`|
| Module: Plotting |[Samtools](https://www.htslib.org/), Gnuplot, matplotlib, tectonic, texlive-core, R and R libraries -- tidyverse ComplexHeatmap, svglite, scales, RColorBrewer|`envs/plot-bamstats.yaml`,`envs/r.yaml`|
</details>

## Steps to take

  * Download the code from this repository.
  * Gather your starting files (see below) and check that they are in the correct format.
  * Edit the `config/config.yaml` to match your files and desired parameters.
  * Instal Mamba/Conda [Microforge3](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
  * Install the `diversity` enviroment: `mamba env create --file workflow/envs/diversity.yaml`
  * `conda activate diversity`
  * Run the pipeline: `snakemake --cores <n> --sdm conda -p`
    * Snakemake options:  
      * `--cores <n>`: Number of cores you want to use. Mandatory.
      * `--sdm conda`: Specify that you will use the conda environments specified in the rules. Mandatory.
      * `-p`: Print the command-lines run by each job in the stdout. Optional.
      * `--conda-frontend conda`: Use it if you are using conda instead of mamba. Mamba is the default.
      * `--rerun-incomplete`: Use it when a past run of the workflow was aborted and you are repeating the run.
      * `--keep-going`: Use it to avoid stoping the workflow when a job fails (i.e. run everything that can run).
      * `-n`: Dry run.

## Structure of the working directory:    
  * `workflow/` has all the code (rules, scripts and main Snakefile), and the environment files.
  * `config/` has the `config.yaml` provided [here](https://github.com/magwenelab/DiversityPipeline/blob/workflow-style/config.yaml) file that **you must edit** to adapt to your dataset.
  * `results/` will hold all the output.
  * `logs/` will hold the log files of all runs.  
  * Additionaly you need to provide the **starting files** described bellow. It s recommended to put the data files in `data/` and the tables in `config/`.

## Starting files:
  * Metadata CSV table: A comma-separated table with one sample per row. Specify the path to it in `config/config.yaml`. Mandatory columns: `sample` (sample ID used in the FASTQ files), `group` (lineage or group that associates the sample with a reference genome), `strain` (strain name, does not need to be unique of each sample, can be the same as `sample`). If plotting will be activated you need  one metadata column to color your samples, specify the name of this column in the `config/config.yaml`. More columns with free format are allowed. [Example](https://github.com/magwenelab/DiversityPipeline/blob/workflow-style/config/sample_metadata.csv)
  * FASTQ files: Paired end short-read FASTQ files, one forward and one reverse file for each sample. The name of these files should be the one used in the metadata `sample` column (followed by an extension specified in the `config/config.yaml`). Files can be gzip compressed. The FASTQ files for all samples should be in the same directory (e.g., `data/samples/`, specified in the `config.yaml`).
  * Reference genomes: FASTA files for each reference genome. The names of the files must be the ones in the `group` column of the metadata, e.g. `VNI.fasta`. If you will use a main reference to annotate the reference genomes provide the FASTA and GFF for the main reference. Otherwise provide the GFF file for the reference genomes. Specify the path for all, and the name of the main reference files in the `config/config.yaml`.

  * `config/loci.csv`: CSV with first column with the gene IDs and second column with the name of the locus/pathway the gene belongs to if you want genes to be plotted to coverage and MAPQ plots.  No column names.[Example](https://github.com/magwenelab/DiversityPipeline/blob/workflow-style/config/loci.csv)
  * `config/chromosome_names.csv`: CSV with colums group, chromosome ID (the sequence ID in the FASTA and GFF of the references), chromosome name (typically a number). No column names. [Example](https://github.com/magwenelab/DiversityPipeline/blob/workflow-style/config/chromosome_names.csv).
  * `config/features.txt`: List of feature names to lift over, more info [here](https://github.com/agshumate/Liftoff?tab=readme-ov-file#feature-types). This file is provided [here](https://github.com/magwenelab/DiversityPipeline/blob/workflow-style/config/features.txt), but you can edit it if you want.


