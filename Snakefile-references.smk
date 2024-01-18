configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
LINS=list(set(samplefile["group"]))

REFDIR = str(config["reference_directory"])
REF_GFF = REFDIR + str(config["reference_gff"])
REF_NAME = str(config["reference_gff"]).replace('.gff', '')

rule all:
    input:
        expand(REFDIR + "{lineage}.gff",lineage=LINS),
        REF_GFF + ".tsv",
        REFDIR + "references_unmapped.svg",

rule ref_gff2tsv:
    input:
        REF_GFF
    output:
        REF_GFF + ".tsv"
    conda:
        "envs/agat.yaml"
    params: 
        REF_NAME    
    log: 
        "logs/references/ref_gff2tsv_agat.log"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input} -o {output} &> {log}"
        " && "
        "rm {params}.agat.log || true" 

rule ref2ref_liftoff:
    input:
        target_refs = REFDIR + "{lineage}.fasta",
        fasta = REFDIR + str(config["reference_fasta"]),
        gff = REF_GFF,
        features = "files/features.txt"
    output:
        lin_gff = REFDIR + "{lineage}.gff",
        unmapped = REFDIR + "{lineage}_unmapped_features.txt"
    threads: config["threads_liftoff"] 
    conda:
        "envs/liftoff.yaml"
    log:
        "logs/references/{lineage}_ref_liftoff.log"
    params:
        refdir = REFDIR    
    shell:
        "liftoff "
        "-g {input.gff} "
        "-polish "
        "-f {input.features} "
        "-o {params.refdir}/{wildcards.lineage}_liftoff.gff "
        "-dir {params.refdir}/{wildcards.lineage}_intermediate_files "
        "-u {output.unmapped} "
        "-p {threads} "
        "{input.target_refs} {input.fasta} "
        "&> {log} "
        "&& "
        "mv {params.refdir}/{wildcards.lineage}_liftoff.gff_polished {output.lin_gff} "
        
rule unmapped_count_plot:
    input:
        REF_GFF + ".tsv",
        config["sample_file"],
        expand(REFDIR + "{lineage}_unmapped_features.txt", lineage=LINS)        
    output:
        REFDIR + "references_unmapped_count.tsv",
        REFDIR + "references_unmapped.svg"
    conda:
        "envs/r.yaml"
    log:
        "logs/references/unmapped_count_plot.log"
    script:
        "scripts/count_reference_unmapped.R"