# Generate a TSV version of the main reference annotation
rule main_gff2tsv: 
    input:
        MAIN_GFF
    output:
        os.path.join(REFDIR, MAIN_NAME + ".tsv")
    conda:
        "../envs/agat.yaml"
    params: 
        ref_name = MAIN_NAME    
    log: 
        "logs/references/ref_gff2tsv_agat.log"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "-gff {input} "
        "-o {output} "
        "&> {log} "
        "&& "
        "rm {params.ref_name}.agat.log || true" 

# Lift over annotation of the main reference to the reference genomes
rule ref2ref_liftoff:
    input:
        target_refs = os.path.join(REFDIR, "{lineage}.fasta"),
        fasta = MAIN_FASTA,
        gff = MAIN_GFF,
        features = FEATURE_FILE
    output:
        target_gff = os.path.join(REFDIR, "{lineage}.gff"),
        unmapped = os.path.join(REFDIR, "{lineage}_unmapped_features.txt")
    threads: 
        config["annotate_references"]["liftoff"]["threads"]
    conda:
        "../envs/liftoff.yaml"
    log:
        "logs/references/{lineage}_ref_liftoff.log"
    params:
        refdir = REFDIR,
        extra = config["annotate_references"]["liftoff"]["extra"]
    shell:
        "liftoff "
        "-g {input.gff} "
        "-polish "
        "-f {input.features} "
        "-o {params.refdir}/{wildcards.lineage}_liftoff.gff "
        "-dir {params.refdir}/{wildcards.lineage}_intermediate_files "
        "-u {output.unmapped} "
        "-p {threads} "
        "{params.extra} "
        "{input.target_refs} {input.fasta} "
        "&> {log} "
        "&& "
        "mv {params.refdir}/{wildcards.lineage}_liftoff.gff_polished {output.target_gff} "
