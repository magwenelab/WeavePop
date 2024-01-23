import os.path
import glob

configfile: "config-annotate-references.yaml"

REF_FASTA = config["main_fasta"]
REF_GFF = config["main_gff"]
REF_NAME, _ = os.path.splitext(os.path.basename(REF_GFF))

SUBREF_DIR = config["subref_directory"]
SUBREF_FASTAS = set(glob.glob(f"{SUBREF_DIR}/*.fasta") + glob.glob(f"{SUBREF_DIR}/*.fa"))
SUBREF_NAMES = [os.path.splitext(os.path.basename(f))[0] for f in SUBREF_FASTAS]

FEATURE_FILE = config["feature_file"]
if not FEATURE_FILE:
    FEATURE_FILE = "features.txt"


rule all:
    input:
        expand(os.path.join(SUBREF_DIR, "{target}.gff"), target=SUBREF_NAMES),
        os.path.join(SUBREF_DIR, REF_GFF + ".tsv"),
        FEATURE_FILE

rule ref_gff2tsv:
    input:
        REF_GFF
    output:
        os.path.join(SUBREF_DIR, REF_GFF + ".tsv")
    conda:
        "envs/agat.yaml"
    params: 
        REF_NAME = REF_NAME    
    log: 
        "logs/references/ref_gff2tsv_agat.log"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input} -o {output} &> {log}"
        " && "
        "rm {params.REF_NAME}.agat.log || true" 


rule ref2ref_liftoff:
    input:
        target_fasta = os.path.join(SUBREF_DIR, "{target}.fasta"),
        ref_fasta = REF_FASTA,
        ref_gff = REF_GFF,
        features = FEATURE_FILE
    output:
        target_gff = os.path.join(SUBREF_DIR, "{target}.gff"),
        target_unmapped = os.path.join(SUBREF_DIR, "{target}_unmapped_features.txt")
    threads: config["threads_liftoff"] 
    conda:
        "envs/liftoff.yaml"
    log:
        "logs/references/{target}_ref_liftoff.log"
    params:
        SUBREF_DIR = SUBREF_DIR
    shell:
        "liftoff "
        "-g {input.ref_gff} "
        "-polish "
        "-f {input.features} "
        "-o {params.SUBREF_DIR}/{wildcards.target}_liftoff.gff "
        "-dir {params.SUBREF_DIR}/{wildcards.target}_intermediate_files "
        "-u {output.target_unmapped} "
        "-p {threads} "
        "{input.target_fasta} {input.ref_fasta} "
        "&> {log} "
        "&& "
        "mv {params.SUBREF_DIR}/{wildcards.target}_liftoff.gff_polished {output.target_gff} "

# rule for generating features.txt if it doesn't exist        
rule generate_features:
    input:
        ref_gff = REF_GFF
    output:
        FEATURE_FILE
    shell:
        "grep -v '^#' {input.ref_gff} | "
        "awk -F'\t' '{{print $3}}' | "
        "sort | "
        "uniq > {output}"