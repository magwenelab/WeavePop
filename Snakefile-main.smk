configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))
ref_table = (pd.read_csv(config["sample_reference_file"], sep=","))
ref_table.set_index('sample', inplace=True)
REFDIR = str(config["reference_directory"])
REF_GFF = REFDIR + str(config["reference_gff"])

protlist=(pd.read_csv("files/protein_list.txt", sep=",", header = None, names = ['protein']))
proteins=list(protlist["protein"])

rule all:
    input:
        expand("analysis/{sample}/snps.consensus.fa",sample=samples),
        expand("analysis/{sample}/snps.bam",sample=samples),
        expand("analysis/{sample}/lifted.gff_polished", sample=samples),
        expand("analysis/{sample}/predicted_cds.fa",sample=samples),
        expand("analysis/{sample}/predicted_proteins.fa",sample=samples),
        "results/proteins.done",
        "results/cds.done",
        "results/unmapped.svg"

rule samples_list:
    output: 
        "files/samples.txt"
    run:
        sample = samplefile["sample"]
        sample.to_csv("files/samples.txt", index = False, header = False)

rule snippy:
    input:
        "fastq_combined/{sample}" + config["fastq_suffix1"],
        "fastq_combined/{sample}" + config["fastq_suffix2"]
    params:
        ref = lambda wildcards: (REFDIR + ref_table.loc[wildcards.sample, 'refgenome']),
        file1 = lambda wildcards: ref_table.loc[wildcards.sample, 'file1'],
        file2 = lambda wildcards: ref_table.loc[wildcards.sample, 'file2'] 
    output:
        "analysis/{sample}/snps.consensus.fa",
        "analysis/{sample}/snps.bam"
    threads: config["threads_snippy"]
    log:
        "logs/snippy/{sample}.log" 
    shell:
        "snippy --outdir analysis/{wildcards.sample} "
        "--cpus {threads} "
        "--ref {params.ref} "
        "--R1 fastq_combined/{params.file1} "
        "--R2 fastq_combined/{params.file2} "
        "--force &> {log}"

rule liftoff:
    input:
        target = "analysis/{sample}/snps.consensus.fa",
        features = "files/features.txt"
    params:
        refgff = lambda wildcards:(REFDIR + ref_table.loc[wildcards.sample, 'group'] + "_liftoff.gff_polished"),
        refgenome = lambda wildcards:(REFDIR + ref_table.loc[wildcards.sample, 'refgenome'])
    output:
        "analysis/{sample}/lifted.gff",        
        "analysis/{sample}/lifted.gff_polished",
        "analysis/{sample}/unmapped_features.txt"
    threads: config["threads_liftoff"]
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
        "rm lifted.agat.log"
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
        expand("analysis/{sample}/predicted_proteins.fa",sample=samples),
        expand("analysis/{sample}/predicted_proteins.fa.fai",sample=samples)
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
        expand("analysis/{sample}/predicted_cds.fa",sample=samples),
        expand("analysis/{sample}/predicted_cds.fa.fai",sample=samples)
    output:
        "results/cds.done"
    log:
        "logs/cds/cds.log"  
    script:
        "scripts/by_cds.sh"

rule unmapped_count_plot:
    input:
        REF_GFF + ".tsv",
        config["sample_file"],
        expand("analysis/{sample}/unmapped_features.txt", sample=samples)        
    output:
        "results/unmapped_count.tsv",
        "results/unmapped.svg"
    log:
        "logs/liftoff/unmapped_count_plot.log"
    script:
        "scripts/count_sample_unmapped.R"
