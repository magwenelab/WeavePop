
# Reference-based mapping, variant calling and assembly / Annotation liftover / Mapping quality and coverage

## Broad description


## Requirements

The environment from which the workflow must mb run has the following software and you can install it with: `mamba env create --file envs/diversity.yml`
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

## Structure of the working directory:    
  * `workflow/` has all the code (rules, scripts and main Snakefile), and the environment files.
  * `config/` has the `config.yaml` provided [here](https://github.com/magwenelab/DiversityPipeline/blob/workflow-style/config.yaml) file that **you must edit** to adapt to your dataset.
  * `results/` will hold all the output.
  * `logs/` will hold the log files of all runs.  
  * Additionaly you need to provide the following **starting files**:
    * Metadata CSV table: A CSV table with one sample per row. Specify the path to it in `config/config.yaml`. Mandatory columns: `sample` (sample ID used in the FASTQ files), `group` (lineage or group that associates the sample with a reference genome), `strain` (strain name, does not be to be unique of each sample). If plotting will be activate you need at least one metadata column to color your samples, specify the name of this column in the `config/config.yaml`. More columns with free format are allowed.
    * FASTQ files: Paired end short-read FASTQ files, one forward and one reverse file for each sample. The name of these files should be the one used in the `sample_metadata.csv` `sample` column (followed by an extension specified in the `config/config.yaml`). Files can be gzip compressed. The FASTQ files for all samples should be in the same directory (e.g., `data/samples/`, specified in the `config.yaml`).
    * Reference genomes: FASTA files for each reference genome.The names of the files must be the ones in the `group` column of the `sample_metadata.csv`, e.g. `VNI.fasta`. If you will use a main reference to annotate the reference genomes provide the FASTA and GFF for the main reference. Otherwise provide the GFF file for the reference genomes. Specify the path and name of the main reference files in the `config/config.yaml`.

    * Optional: `config/loci.csv`, CSV with first column with the gene IDs and second column with the name of the locus/pathway the gene belongs to if you want genes to be plotted to coverage and MAPQ plots  No column names.
    * `config/chromosome_names.csv` with columns (without column names): group, chromosome ID (the sequence ID in the Fasta and GFF of the references), chromosome name (typically a number). If your genomes are Complete Genomes from NCBI use `bash get-chromosome_names.sh` to get this file.
    * `config/features.txt` list of feature names to lift over. This file is provided in this repository.
    * `references/` directory with:
      * Fasta files to use as reference genomes. The names of the files must be the ones in the "group" column of the `files/sample_metadata.csv`, e.g. `VNI.fasta`
      * Optional: Fasta and GFF files of main reference (one with available annotation with the desired gene IDs). If your main reference is one of your reference genomes, duplicate the genome files and give another name to the ones that will be used as main.
      * If main reference is not provided, GFF files of reference genomes are needed (with the same name as the fastas).
        * If your genomes have a mitochondrial chromosome you can run `bash get-removed-chrom.sh path-to-fasta path-to-gff seq_id` to remove it, in an environment with Seqkit available (diversity).

  See example files [here](https://github.com/magwenelab/DiversityPipeline/tree/workflow-style/config).

