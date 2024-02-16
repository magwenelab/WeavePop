# Generate softlinks of main reference
rule main_links:
    input:
        gff = MAIN_GFF,
        fasta = MAIN_FASTA
    output:
        gff = REFDIR / "{lineage}" / str(MAIN_NAME + ".gff"),
        fasta = REFDIR / "{lineage}" / str(MAIN_NAME + ".fasta")
    log:
        "logs/references/main_liks_{lineage}.log"
    shell:
        "ln -s -r {input.gff} {output.gff} &> {log}"
        "&& "
        "ln -s -r {input.fasta} {output.fasta} &>> {log}"

# Generate a TSV version of the main reference annotation
rule main_gff2tsv: 
    input:
        MAIN_GFF
    output:
        REFDIR / str(MAIN_NAME + ".tsv")
    conda:
        "../envs/agat.yaml"
    params: 
        ref_name = MAIN_NAME    
    log: 
        "logs/references/main_gff2tsv.log"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "-gff {input} "
        "-o {output} "
        "&> {log} "
        "&& "
        "rm {params.ref_name}.agat.log &>> {log} || true" 

# Lift over annotation of the main reference to the reference genomes
rule ref2ref_liftoff:
    input:
        target_refs = REFDIR / "{lineage}" / "{lineage}.fasta",
        fasta = rules.main_links.output.fasta,
        gff = rules.main_links.output.gff,
        features = FEATURE_FILE
    output:
        target_gff = REFDIR / "{lineage}" / "{lineage}.gff",
        unmapped = REFDIR / "{lineage}" / "unmapped_features.txt"
    threads: 
        config["annotate_references"]["liftoff"]["threads"]
    conda:
        "../envs/liftoff.yaml"
    log:
        "logs/references/ref2ref_liftoff_{lineage}.log"
    params:
        refdir = REFDIR,
        extra = config["annotate_references"]["liftoff"]["extra"]
    shell:
        "liftoff "
        "-g {input.gff} "
        "-o {params.refdir}/{wildcards.lineage}/liftoff.gff "
        "-dir {params.refdir}/{wildcards.lineage}/intermediate_files "
        "-u {output.unmapped} "
        "-p {threads} "
        "-f {input.features} "
        "-polish "
        "{params.extra} "
        "{input.target_refs} {input.fasta} "
        "&> {log} "
        "&& "
        "mv {params.refdir}/{wildcards.lineage}/liftoff.gff_polished {output.target_gff} &>> {log}"