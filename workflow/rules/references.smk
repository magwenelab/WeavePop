# =================================================================================================
#   Per main reference | Standardize GFF format and convert to TSV
# =================================================================================================

# Run AGAT to add and modify tags and convert to TSV 
rule agat_fix_gff:
    input: 
        gff= MAIN_GFF,
        config = rules.agat_config.output
    output:
        fixed_ID = os.path.join(INT_REFS_DIR, f"{MAIN_NAME}_fixed_ID.gff"),
        fixed_locus = os.path.join(INT_REFS_DIR, f"{MAIN_NAME}_fixed_locus.gff"),
        fixed_description = os.path.join(INT_REFS_DIR, f"{MAIN_NAME}_fixed_description.gff"),
        tsv = os.path.join(INT_REFS_DIR, f"{MAIN_NAME}_fixed.tsv")
    params:
        name = MAIN_NAME
    log:
        "logs/references/agat_fix_gff.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/agat.yaml"
    shell:
        """
        agat_convert_sp_gxf2gxf.pl -g {input.gff} -o {output.fixed_ID} -c {input.config} &> {log} && 
        agat_sq_add_locus_tag.pl --gff {output.fixed_ID} --li ID -o {output.fixed_locus} -c {input.config} &>> {log} && 
        agat_sp_manage_attributes.pl --gff {output.fixed_locus} --tag product/description -o {output.fixed_description}  -c {input.config} &>> {log} && 
        agat_convert_sp_gff2tsv.pl --gff {output.fixed_description} -o {output.tsv} -c {input.config} &>> {log} 
        """

# Recreate IDs
rule fix_gff_tsv:
    input:
        tsv = rules.agat_fix_gff.output.tsv
    output:
        gff = os.path.join(INT_REFS_DIR, f"{MAIN_NAME}.gff"),
        tsv = os.path.join(INT_REFS_DIR, f"{MAIN_NAME}.tsv")
    log:
        "logs/references/fix_gff_tsv.log"
    conda:
        "../envs/snakemake.yaml"
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
        fasta = os.path.join(INT_REFS_DIR, "{lineage}", f"{MAIN_NAME}.fasta"),
        gff = os.path.join(INT_REFS_DIR, "{lineage}", f"{MAIN_NAME}.gff")
    log:
        "logs/references/main_liks_{lineage}.log"
    conda:
        "../envs/shell.yaml"
    shell:
        "ln -s -r {input.fasta} {output.fasta} &>> {log}"
        "&& "
        "ln -s -r {input.gff} {output.gff} &>> {log}"

# ==================================================================================================
#   Per lineage | Annotate reference genomes
# ==================================================================================================

# Run Lifotff to annotate reference genomes with main reference
rule ref2ref_liftoff:
    input:
        target_refs = INT_REFS_DIR / "{lineage}" / "{lineage}.fasta",
        fasta = rules.main_links.output.fasta,
        gff = rules.main_links.output.gff
    output:
        target_gff = REFS_DIR / "{lineage}.gff",
        unmapped = INT_REFS_DIR / "{lineage}" / "unmapped_features.txt",
        intermediate = directory(INT_REFS_DIR / "{lineage}" / "intermediate_liftoff"),
        fai_main = os.path.join(INT_REFS_DIR, "{lineage}", f"{MAIN_NAME}.fasta.fai"),
        fai = INT_REFS_DIR / "{lineage}" / "{lineage}.fasta.fai",
        mmi = INT_REFS_DIR / "{lineage}" / "{lineage}.fasta.mmi"
    params:
        refdir = INT_REFS_DIR,
        extra = config["annotate_references"]["liftoff"]["extra"]
    log:
        "logs/references/ref2ref_liftoff_{lineage}.log"
    threads: 
        config["annotate_references"]["liftoff"]["threads"]
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/liftoff.yaml"
    shell:
        "liftoff "
        "-g {input.gff} "
        "-o {params.refdir}/{wildcards.lineage}/liftoff.gff "
        "-dir {output.intermediate} "
        "-u {output.unmapped} "
        "-p {threads} "
        "-polish "
        "{params.extra} "
        "{input.target_refs} {input.fasta} "
        "&> {log} "
        "&& "
        "mv {params.refdir}/{wildcards.lineage}/liftoff.gff_polished {output.target_gff} &>> {log}"

# Convert GFF file to TSV format
rule gff2tsv:
    input:
        target = rules.ref2ref_liftoff.output.target_gff,
        config = rules.agat_config.output
    output:
        INT_REFS_DIR / "{lineage}" / "{lineage}.gff.tsv"
    log:
        "logs/references/gff2tsv_{lineage}.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input.target} -c {input.config} -o {output} &> {log} "
