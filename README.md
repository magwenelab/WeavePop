
# Variant-calling / reference-based assembly / annotation liftover / quality and coverage analysis pipeline

## Broad description


## Requirements

* Miniconda3
* Python
* Python modules -- Pandas, Scipy
* Xonsh -- https://xon.sh/
* R with tidyverse meta-package
* Snakemake -- https://snakemake.github.io/
  * Graphviz -- https://graphviz.org/ (optional, to see Snakemake DAG in a graph) 
* NCBI Entrez Utilities (E-utilities) command line tools -- https://www.ncbi.nlm.nih.gov/books/NBK25501/
* NCBI SRA Tools -- https://github.com/ncbi/sra-tools
* Seqkit -- https://bioinf.shenwei.me/seqkit/
* Snippy -- https://github.com/tseemann/snippy
* Liftoff -- https://github.com/agshumate/Liftoff
* AGAT -- https://github.com/NBISweden/AGAT
* Mosdepth -- https://github.com/brentp/mosdepth
* Samtools -- https://www.htslib.org/
  
### Installations  

Environment installation files are in `envs/`
<details>
<summary>Install crypto_div environment -- everything runs in this environment </summary>
 
conda env create --file envs/crypto_div.yml

</details>


<details>
<summary>Install depth environment </summary>
 
conda env create --file envs/depth.yml
</details>

<details>
<summary> Install AGAT in its own environment -- when used in Snakemake the workflow uses a separate installation of the environment, but to use it independently install it with this instructions </summary>

Run this lines one by one:
~~~
 conda create -n agat
 conda activate agat
 conda install perl-bioperl perl-clone perl-graph perl-lwp-simple perl-carp perl-sort-naturally perl-file-share perl-file-sharedir-install perl-moose perl-yaml perl-lwp-protocol-https -c bioconda
 conda install r-base
 conda install perl-statistics-r -c bioconda
 cpan install bioperl List::MoreUtils Term::ProgressBar
 git clone https://github.com/NBISweden/AGAT.git
 perl Makefile.PL 
 make
 make test
 make install
 conda deactivate
 ~~~

</details>


## Overview

### Starting files: 
  * `files/sample_metadata.csv` with columns: strain, sample (the one in the fastq file names), group (lineage or group to associate to a reference genome), more-optional-metadata-fields
  * `files/lineage_references.csv` with columns: group, file (file name of reference genome assembly), strain, more-optional-metadata-fields (like genbank accession and bioproject)
   * Lists of genes of loci of interest:
    * `files/centromere.txt` (with IDs of genes in main reference GFF)
  * `files/chromosome_names.csv` with columns: group, chromosome ID (the sequence ID in the Fasta and GFF of the references), chromosome name (typically a number). Without column names. Use `get-chromosome_names.sh` if your genomes are Complete Genomes from NCBI to get this file.
  * `references/` directory with:
    * Fasta files to use as reference genomes.
    * Fasta file of main reference (the one with available annotation with the desired gene IDs)
    * GFF file of main reference
 

### Structure of repository:
  * The working directory has the scripts and Snakefiles to run.  
  * `files/` has some of the starting files and files created by the pipeline.
  * `scripts/` has the scripts used by the Snakefiles, not by the user directly.  
  * `analyses/` has one directory per sample, all the resulting files of the analyses performed per sample are there with a generic name.  
  * `results/` has the resulting files of the analyses that consider all the samples.  

### Scripts to be run in this order:

#### Module 1 (Optional): To download all fastqs of a BioProject
1. Get files: `xonsh scripts/get-seqdata-of-bioproject.xsh -p PRJNA685103`  
2. Combine fastqs of the same sample, rename with sample ID and compress:
   `parallel xonsh scripts/fastq-combiner.xsh {} files/read_pair_table.csv fastqs/ fastq_combined/ :::: files/samples.txt`
   
#### Module 1: Annotate references according to main reference
`Snakefile-references.smk` -- is a Snakefile to lift over annotations from `FungiDB-53_CneoformansH99_PMM.gff` into the four lineage genomes (`{lineage}_{GenBank Accession}.fasta`).  
    * It currently works with:  
  ` snakemake --snakefile Snakefile-references.smk --cores 1 --use-conda --conda-frontend conda -p`:  
      ⚠️ `--cores 1` is because there is a problem if Liftoff runs in parallel because the different jobs try to create `FungiDB-65_CneoformansH99.gff_db` at the same time and that is not cool.    
      ⚠️ `--conda-frontend conda` because it cannot use mamba, which is the default.  
      ⏰ Pending: Merge into main workflow.
    * Output:  
  
      *  `references/reference_genes.tsv`
      *  `references/{lineage}_liftoff.gff_polished`
      *  `references/{lineage}_liftoff.gff_polished.tsv`
      *  `references/{lineage}_predicted_proteins.fa`
      *  `references/{lineage}_predicted_cds.fa`
      *  `results/protein_list.txt`
      *  `references/references_unmapped_features.csv`
      *  `references/references_unmapped_count.csv`
      *  `references/references_unmapped.png`
      * And more intermediate and extra files

#### Module 2: Main analyses
`Snakefile-main.smk`-- is the Snakefile to run the pipeline, it uses the `config.yaml` file.  
It runs the script `scripts/fastq-combiner.xsh` for each sample in `files/read_pair_table.csv`. This concatenates all `_1.fastq` of one sample into only one file named `{SRS-accession}_1.fq.gz` and compresses it and does the same for `_2.fastq`.  
It runs **snippy**, **liftoff** and **agat** for each sample, it **extracts sequences** (cds and protein) of each sample and **concatenates** them by cds and by protein.

    * Output:  
    
      * `fastq_combined/{SRS-accession}_1.fq.gz` and `fastq_combined/{SRS-accession}_2.fq.gz`.
      * `analyses/{sample}/snps.consensus.fa` and extra assembly files    
      * `analyses/{sample}/snps.bam` and extra alignment files  
      * `analyses/{sample}/snps.vcf` and extra variant calling files  
      * `analyses/{sample}/lifted.gff_polished` and extra annotation files  
      * `analyses/{sample}/predicted_cds.fa`  and extra index files
      * `analyses/{sample}/predicted_proteins.fa`  and extra index files
      * `results/cds/{protein}.fa`
      * `results/proteins/{protein}.fa`

#### Module 3: Quality and depth analyses
`Snakefile-depth-quality.smk`: Generates **quality and coverage** plots.  
   * Output:  
  
     * `results/mapped_reads.svg` and `results/mapping_stats.txt` plot and table with fraction of mapping reads per sample.  
     * `analyses/{sample}/snps.bam.stats` and `analyses/{sample}/bamstats/` directory with `plot-bamstats` resulting plots.  
     * `analyses/{sample}/coverage.svg` depth of coverage along chromosome plot (qith all and only good quality mappings and location of interesting loci).  
     * `analyses/{sample}/coverage_stats.svg` mean and median coverage pero chromosome and global.  
     * `analyses/{sample}/cov_distribution_.svg` ditribution of coverage values plot.  
     * `analyses/{sample}/mapq.svg` mapping quality along chromosome plot.  
     * `analyses/{sample}/mapq_distribution.svg` distribution of maping quality values plot.    
     * `analyses/{sample}/annotation.gff` GFF file with complete annotation plus average MAPQ and coverage of windows in which the features are located.    
 
