configfile: "config.yaml"

import pandas as pd

lin_ref_table = (pd.read_csv(config["lineage_reference_file"], sep=","))
LINS=list(lin_ref_table["group"])
lin_ref_table.set_index('group', inplace=True)

REFDIR = str(config["reference_directory"])
REF_GFF = REFDIR + str(config["reference_gff"])

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
    log: 
        "logs/references/ref_gff2tsv_agat.log"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input} -o {output} &> {log}"

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
        config["lineage_reference_file"],
        expand(REFDIR + "{lineage}_unmapped_features.txt", lineage=LINS)        
    output:
        REFDIR + "references_unmapped_count.tsv",
        REFDIR + "references_unmapped.svg"
    log:
        "logs/references/unmapped_count_plot.log"
    script:
        "scripts/count_reference_unmapped.R"