configfile: "config.yaml"

import pandas as pd

lin_ref_table = (pd.read_csv(config["lineage_reference_file"], sep=","))
LINS=list(lin_ref_table["group"])
lin_ref_table.set_index('group', inplace=True)

REFDIR = str(config["reference_directory"])
REF_GFF = REFDIR + str(config["reference_gff"])

rule all:
    input:
        expand(REFDIR + "{lineage}_liftoff.gff_polished",lineage=LINS),
        expand(REFDIR + "{lineage}_liftoff.gff_polished.tsv",lineage=LINS),
        expand(REFDIR + "{lineage}_predicted_proteins.fa",lineage=LINS),
        expand(REFDIR + "{lineage}_predicted_cds.fa",lineage=LINS),
        "files/protein_list.txt",
        REFDIR + "references_unmapped.svg",
        config["locitsv"],
        config["sample_reference_file"]

rule features:
    input:
        REF_GFF
    output:
        REFDIR + "features.txt"
    log:
        "logs/references/features.log"
    shell:
        "grep -v '#' {input} | cut -f 3 | sort | uniq > {output} 2> {log}"
rule ref2ref_liftoff:
    input:
        target_refs = REFDIR + "{lineage}.fasta",
        fasta = REFDIR + str(config["reference_fasta"]),
        gff = REF_GFF,
        features = REFDIR + "features.txt"
    output:
        lin_gff = REFDIR + "{lineage}_liftoff.gff_polished",
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
rule gff2tsv:
    input:
        REFDIR + "{lineage}_liftoff.gff_polished"
    output:
        REFDIR + "{lineage}_liftoff.gff_polished.tsv",
        temp("{lineage}_liftoff.agat.log")
    conda:
        "envs/agat.yaml"
    log:
        "logs/references/{lineage}_gff2tsv.log"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input} -o {output} "
        "&> {log} "

rule loci_interest:
    input:
        expand(REFDIR + "{lineage}_liftoff.gff_polished.tsv", lineage=LINS)
    output:
        config["locitsv"]
    params:
        loci=config["loci"]
    log: 
        "logs/references/loci.log"
    shell:
        "xonsh scripts/loci.xsh {params.loci} -o {output} {input} &> {log}"

rule ref2ref_agat:
    input: 
        lin_liftoff = REFDIR + "{lineage}_liftoff.gff_polished",
        lin_fasta = REFDIR + "{lineage}.fasta"
    output:
        temp("{lineage}_liftoff.agat.log"),
        cds = REFDIR + "{lineage}_predicted_cds.fa",
        prots = REFDIR + "{lineage}_predicted_proteins.fa"
    conda:
        "envs/agat.yaml"
    log:
        cds = "logs/references/{lineage}_ref_agat_cds.log",   
        prots = "logs/references/{lineage}_ref_agat_prots.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.lin_liftoff} "
        "-f {input.lin_fasta} "
        "-o {output.cds} "
        "&> {log.cds} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input.lin_liftoff} "
        "-f {input.lin_fasta} "
        "-o {output.prots} "
        "-p &> {log.prots}"  

rule protein_list:
    input:
        fasta = REFDIR + "{lineage}_predicted_proteins.fa"
    output:
        list = REFDIR + "{lineage}_protein_list.txt"
    log:
        "logs/references/{lineage}_protein_list.log"
    shell:
        "seqkit seq -n -i {input.fasta} > {output.list} 2> {log}"

rule cat_lists:
    input: 
        expand(REFDIR + "{lineage}_protein_list.txt", lineage=LINS)
    output:
        "files/protein_list.txt"
    log:
        "logs/references/cat_list.log"
    shell:
        "cat {input} | sort | uniq > {output} 2> {log}"
rule ref_gff2tsv:
    input:
        REF_GFF
    output:
        complete = REF_GFF + ".tsv"
    conda:
        "envs/agat.yaml"
    log: 
        "logs/references/ref_gff2tsv_agat.log"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input} -o {output.complete} &> {log}"

rule unmapped_count_plot:
    input:
        REF_GFF + ".tsv",
        config["lineage_reference_file"],
        expand(REFDIR + "{lineage}_unmapped_features.txt", lineage=LINS)        
    output:
        REFDIR + "references_unmapped_count.txt",
        REFDIR + "references_unmapped.svg"
    log:
        "logs/references/unmapped_count_plot.log"
    script:
        "scripts/count_reference_unmapped.R"

rule reference_table:
    input:
        s = config["sample_file"],
        r = config["lineage_reference_file"]  
    output:
        config["sample_reference_file"]
    log:
        "logs/references/reftable.log"
    shell:
        "xonsh scripts/reference_table.xsh -s {input.s} -r {input.r} -o {output} &> {log}"