rule agat_fix_gff:
    input: 
        gff= MAIN_GFF,
        config = rules.agat_config.output
    output:
        fixed_ID = temp(REFDIR / str(MAIN_NAME + "_fixed_ID.gff")),
        fixed_locus = temp(REFDIR / str(MAIN_NAME + "_fixed_locus.gff")),
        fixed_description = temp(REFDIR / str(MAIN_NAME + "_fixed_description.gff")),
        tsv = temp(REFDIR / str(MAIN_NAME + "_fixed.tsv"))
    conda:
        "../envs/agat.yaml"
    params:
        name = MAIN_NAME
    log:
        "logs/references/agat_fix_gff.log"
    shell:
        """
        agat_convert_sp_gxf2gxf.pl -g {input.gff} -o {output.fixed_ID} -c {input.config} &> {log} && 
        agat_sq_add_locus_tag.pl --gff {output.fixed_ID} --li ID -o {output.fixed_locus} -c {input.config} &>> {log} && 
        agat_sp_manage_attributes.pl --gff {output.fixed_locus} --tag product/description -o {output.fixed_description}  -c {input.config} &>> {log} && 
        agat_convert_sp_gff2tsv.pl --gff {output.fixed_description} -o {output.tsv} -c {input.config} &>> {log} 
        """

rule fix_gff_tsv:
    input:
        tsv = rules.agat_fix_gff.output.tsv
    output:
        gff = REFDIR / str(MAIN_NAME + ".gff"),
        tsv = REFDIR / str(MAIN_NAME + ".tsv")
    log:
        "logs/references/fix_gff_tsv.log"
    shell:
        """
        python workflow/scripts/fix_gff.py -i {input.tsv} -og {output.gff} -ot {output.tsv} &> {log}
        """

# Generate softlinks of main reference
rule main_links:
    input:
        fasta = MAIN_FASTA,
        gff = rules.fix_gff_tsv.output.gff
    output:
        fasta = REFDIR / "{lineage}" / str(MAIN_NAME + ".fasta"),
        gff = REFDIR / "{lineage}" / str(MAIN_NAME + ".gff")
    log:
        "logs/references/main_liks_{lineage}.log"
    shell:
        "ln -s -r {input.fasta} {output.fasta} &>> {log}"
        "&& "
        "ln -s -r {input.gff} {output.gff} &>> {log}"

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

# Convert GFF file to TSV format
rule gff2tsv:
    input:
        rules.ref2ref_liftoff.output.target_gff
    output:
        REFDIR / "{lineage}" / "{lineage}.gff.tsv"
    conda:
        "../envs/agat.yaml"
    log:
        "logs/references/gff2tsv_{lineage}.log"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input} -o {output} "
        "&> {log} && "
        "rm {wildcards.lineage}.agat.log &>> {log} || true"