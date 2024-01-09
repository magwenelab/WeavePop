configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
SAMPLES=list(set(samplefile["sample"]))
REFDIR = str(config["reference_directory"])
lin_ref_table = (pd.read_csv(config["lineage_reference_file"], sep=","))
LINS=list(lin_ref_table["group"])

rule all:
    input:
        expand(REFDIR + "{lineage}.gff.tsv", lineage=LINS),
        config["locitsv"],
        expand(REFDIR + "{lineage}_predicted_cds.fa", lineage=LINS),
        expand(REFDIR + "{lineage}_predicted_proteins.fa", lineage=LINS),
        expand(REFDIR + "{lineage}_protein_list.txt", lineage=LINS),
        expand("analysis/{sample}/snps.consensus.fa",sample=SAMPLES),
        expand("analysis/{sample}/snps.bam",sample=SAMPLES),
        expand("analysis/{sample}/lifted.gff_polished", sample=SAMPLES),
        expand("analysis/{sample}/predicted_cds.fa",sample=SAMPLES),
        expand("analysis/{sample}/predicted_proteins.fa",sample=SAMPLES),
        "results/proteins.done",
        "results/cds.done",
        # "results/unmapped.svg"

rule samples_list:
    output: 
        "files/samples.txt"
    run:
        sample = samplefile["sample"]
        sample.to_csv("files/samples.txt", index = False, header = False)

rule reference_table:
    input:
        s = config["sample_file"],
        r = config["lineage_reference_file"]  
    output:
        config["sample_reference_file"]
    params:
        f1 = config["fastq_suffix1"],
        f2= config["fastq_suffix2"] 
    log:
        "logs/reftable.log"
    shell:
        "xonsh scripts/reference_table.xsh -s {input.s} -r {input.r} -o {output} -f1 {params.f1} -f2 {params.f2} &> {log}"

rule gff2tsv:
    input:
        REFDIR + "{lineage}.gff"
    output:
        REFDIR + "{lineage}.gff.tsv"
    conda:
        "envs/agat.yaml"
    log:
        "logs/references/{lineage}_gff2tsv.log"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input} -o {output} "
        "&> {log} && "
        "rm {wildcards.lineage}.agat.log"

rule loci:
    input:
        expand(REFDIR + "{lineage}.gff.tsv", lineage=LINS)
    output:
        config["locitsv"]
    params:
        loci=config["loci"]
    log: 
        "logs/references/loci.log"
    run:
        if config["loci"] == "":
            shell("touch {output}")
        else:
            shell("xonsh scripts/loci.xsh {params.loci} -o {output} {input} &> {log}")

rule ref_agat:
    input: 
        lin_gff = REFDIR + "{lineage}.gff",
        lin_fasta = REFDIR + "{lineage}.fasta"
    output:
        cds = REFDIR + "{lineage}_predicted_cds.fa",
        prots = REFDIR + "{lineage}_predicted_proteins.fa"
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
        "rm {wildcards.lineage}.agat.log"  

rule protein_list:
    input:
        fasta = REFDIR + "{lineage}_predicted_proteins.fa"
    output:
        list = REFDIR + "{lineage}_protein_list.txt"
    conda:
        "envs/agat.yaml"
    log:
        "logs/references/{lineage}_protein_list.log"
    shell:
        "seqkit seq -n -i {input.fasta} 1> {output.list} 2> {log}"

rule cat_lists:
    input: 
        expand(REFDIR + "{lineage}_protein_list.txt", lineage=LINS)
    output:
        "files/protein_list.txt"
    log:
        "logs/references/cat_list.log"
    shell:
        "cat {input} | sort | uniq > {output} 2> {log}"

rule snippy:
    input:
        config["fastq_directory"] + "{sample}" + config["fastq_suffix1"],
        config["fastq_directory"] + "{sample}" + config["fastq_suffix2"],
        config["sample_reference_file"]
    params:
        ref = lambda wildcards: (REFDIR + pd.read_csv(config["sample_reference_file"], sep = ",", index_col=['sample']).loc[wildcards.sample,'refgenome']),
        file1 = lambda wildcards: (pd.read_csv(config["sample_reference_file"], sep = ",", index_col=['sample']).loc[wildcards.sample,'file1']),
        file2 = lambda wildcards: (pd.read_csv(config["sample_reference_file"], sep = ",", index_col=['sample']).loc[wildcards.sample,'file2']),
        fqdir = config["fastq_directory"] 
    output:
        "analysis/{sample}/snps.consensus.fa",
        "analysis/{sample}/snps.bam"
    threads: 
        config["threads_snippy"]
    log:
        "logs/snippy/{sample}.log"
    shell:
        "snippy --outdir analysis/{wildcards.sample} "
        "--cpus {threads} "
        "--ref {params.ref} "
        "--R1 {params.fqdir}{params.file1} "
        "--R2 {params.fqdir}{params.file2} "
        "--force &> {log}"

rule liftoff:
    input:
        config["sample_reference_file"],
        target = "analysis/{sample}/snps.consensus.fa",
        features = "files/features.txt"
    params:
        refgff = lambda wildcards:(REFDIR + pd.read_csv(config["sample_reference_file"], sep = ",", index_col=['sample']).loc[wildcards.sample, 'group'] + ".gff"),
        refgenome = lambda wildcards:(REFDIR + pd.read_csv(config["sample_reference_file"], sep = ",", index_col=['sample']).loc[wildcards.sample, 'refgenome'])
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

# rule unmapped_count_plot:
#     input:
#         REF_GFF + ".tsv",
#         config["sample_file"],
#         expand("analysis/{sample}/unmapped_features.txt", sample=SAMPLES)        
#     output:
#         "results/unmapped_count.tsv",
#         "results/unmapped.svg"
#     log:
#         "logs/liftoff/unmapped_count_plot.log"
#     script:
#         "scripts/count_sample_unmapped.R"
