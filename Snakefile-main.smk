import os.path
import pandas as pd


configfile: "config-main.yaml"

SAMPLE_TABLE = pd.read_csv(config["sample_table"]).set_index("sample", drop=False)
SAMPLES = list(set(SAMPLETBL["sample"]))
LINEAGES=list(set(SAMPLETBL["lineages"]))

REFDIR = config["reference_directory"]
OUTDIR = config["output_directory"]
FASTQDIR = config["fastq_directory"]


rule all:
    input:
        expand(os.path.join(REFDIR, "{lineage}_predicted_cds.fa"), lineage=LINEAGES),
        expand(os.path.join(REFDIR, "{lineage}_predicted_proteins.fa"), lineage=LINEAGES),
        expand(os.path.join(REFDIR, "{lineage}_protein_list.txt"), lineage=LINEAGES),
        expand(os.path.join(OUTDIR, "{sample}/snps.consensus.fa"), sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/snps.bam"), sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/lifted.gff_polished"), sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/predicted_cds.fa"), sample=SAMPLES),
        expand(os.path.join(OUTDIR, "{sample}/predicted_proteins.fa"), sample=SAMPLES),
        os.path.join(OUTDIR, "proteins.done"),
        os.path.join(OUTDIR, "cds.done"),


rule samples_list:
    output: 
        os.path.join(OUTDIR, "samples.txt")
    run:
        sample = SAMPLETBL["sample"]
        sample.to_csv("{output}", index = False, header = False)

# rule reference_table:
#     input:
#         config["sample_file"],
#     output:
#         "files/sample_reference.csv"
#     params:
#         f1 = config["fastq_suffix1"],
#         f2= config["fastq_suffix2"] 
#     log:
#         "logs/reftable.log"
#     shell:
#         "scripts/build_reference_table.py -s {input} -o {output} -f1 {params.f1} -f2 {params.f2} &> {log}"


rule ref_agat:
    input: 
        lin_gff = os.path.join(REFDIR, "{lineage}.gff"),
        lin_fasta = os.path.join(REFDIR, "{lineage}.fasta")
    output:
        cds = os.path.join(REFDIR,"{lineage}_predicted_cds.fa"),
        prots = os.path.join(REFDIR,"{lineage}_predicted_proteins.fa")
    conda:
        "envs/agat.yaml"
    log:
        cds = "logs/references/{lineage}_ref_agat_cds.log",   
        prots = "logs/references/{lineage}_ref_agat_prots.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.lin_gff} "
        "-f {input.lin_fasta} "
        "-o {output.cds} "
        "&> {log.cds} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input.lin_gff} "
        "-f {input.lin_fasta} "
        "-o {output.prots} "
        "-p &> {log.prots} "
        " && "
        "rm {wildcards.lineage}.agat.log || true"  

rule ref_protein_list:
    input:
        fasta = os.path.join(REFDIR, "{lineage}_predicted_proteins.fa")
    output:
        protein_list = os.path.join(REFDIR, "{lineage}_protein_list.txt")
    conda:
        "envs/agat.yaml"
    log:
        "logs/references/{lineage}_protein_list.log"
    shell:
        "seqkit seq -n -i {input.fasta} 1> {output.protein_list} 2> {log}"


rule cat_protein_lists:
    input: 
        expand(os.path.join(REFDIR, "{lineage}_protein_list.txt"), lineage=LINS)
    output:
        "files/protein_list.txt"
    log:
        "logs/references/cat_list.log"
    shell:
        "cat {input} | sort | uniq > {output} 2> {log}"


def fq1_from_sample(wildcards):
  return SAMPLE_TABLE.loc[wildcards.sample, "fastq1"]

def fq2_from_sample(wildcards):
  return SAMPLE_TABLE.loc[wildcards.sample, "fastq2"]

def lineage_from_sample(wildcards):
    return SAMPLE_TABLE.loc[wildcards.sample, "lineage"]

def refgenome_from_sample(wildcards):
    return SAMPLE_TABLE.loc[wildcards.sample, "refgenome"]


rule snippy:
    input:
        fq1 = fq1_from_sample,
        fq2 = fq2_from_sample,
        refgenome = refgenome_from_sample,
        lineage = lineage_from_sample 
    params:
        outdir = OUTDIR
        fqdir = FASTQDIR
    output:
        os.path.join(OUTDIR, "{sample}/snps.consensus.fa"),
        os.path.join(OUTDIR, "{sample}/snps.bam")
    threads: 
        config["threads_snippy"]
    conda:
        "envs/snippy.yaml"
    log:
        "logs/snippy/{sample}.log"
    shell:
        "snippy --outdir {OUTDIR}/{wildcards.sample} "
        "--cpus {threads} "
        "--ref {input.refgenome} "
        "--R1 {FASTQDIR}/{input.fq1} "
        "--R2 {FASTQDIR}/{input.fq2} "
        "--force &> {log}"


rule liftoff:
    input:
        "files/sample_reference.csv",
        target = "analysis/{sample}/snps.consensus.fa",
        features = "files/features.txt"
    params:
        refgff = lambda wildcards:(REFDIR + pd.read_csv("files/sample_reference.csv", sep = ",", index_col=['sample']).loc[wildcards.sample, 'group'] + ".gff"),
        refgenome = lambda wildcards:(REFDIR + pd.read_csv("files/sample_reference.csv", sep = ",", index_col=['sample']).loc[wildcards.sample, 'refgenome'])
    output:
        "analysis/{sample}/lifted.gff",        
        "analysis/{sample}/lifted.gff_polished",
        "analysis/{sample}/unmapped_features.txt"
    threads: config["threads_liftoff"]
    conda:
        "envs/liftoff.yaml"
    log:
        "logs/liftoff/{sample}.log" 
    shell:
        "liftoff "
        "-g {params.refgff} " 
        "-polish "
        "-f {input.features} "
        "-dir analysis/{wildcards.sample}/intermediate_files "
        "-u analysis/{wildcards.sample}/unmapped_features.txt "
        "-o analysis/{wildcards.sample}/lifted.gff "
        "-p {threads} "
        "{input.target} "
        "{params.refgenome} &> {log}"

rule agat:
    input:
        gff = "analysis/{sample}/lifted.gff_polished",
        fa = "analysis/{sample}/snps.consensus.fa"
    output:
        cds = "analysis/{sample}/predicted_cds.fa",
        prots = "analysis/{sample}/predicted_proteins.fa"
    conda:
        "envs/agat.yaml"
    log: 
        cds = "logs/agat/{sample}_cds.log",
        prots = "logs/agat/{sample}_prots.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.cds} "
        "&> {log.cds} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.prots} "
        "-p  &> {log.prots} " 
        " && "
        "rm lifted.agat.log || true"

rule index_proteins:
    input:
        "analysis/{sample}/predicted_proteins.fa"
    output:
        "analysis/{sample}/predicted_proteins.fa.fai"
    conda:
        "envs/agat.yaml"        
    log:
        "logs/faidx/{sample}_proteins.log"    
    shell:
        "seqkit faidx {input} &> {log}"

rule index_cds:
    input:
        "analysis/{sample}/predicted_cds.fa"
    output:
        "analysis/{sample}/predicted_cds.fa.fai"
    conda:
        "envs/agat.yaml"        
    log:
        "logs/faidx/{sample}_cds.log"      
    shell:
        "seqkit faidx {input} &> {log}"

rule by_proteins:
    input:
        "files/protein_list.txt",
        "files/samples.txt",
        expand("analysis/{sample}/predicted_proteins.fa",sample=SAMPLES),
        expand("analysis/{sample}/predicted_proteins.fa.fai",sample=SAMPLES)
    output:
        "results/proteins.done"
    log:
        "logs/proteins/proteins.log"  
    script:
        "scripts/by_proteins.sh"

rule by_cds:
    input:
        "files/protein_list.txt",
        "files/samples.txt",
        expand("analysis/{sample}/predicted_cds.fa",sample=SAMPLES),
        expand("analysis/{sample}/predicted_cds.fa.fai",sample=SAMPLES)
    output:
        "results/cds.done"
    log:
        "logs/cds/cds.log"  
    script:
        "scripts/by_cds.sh"

