# =================================================================================================
#   Main reference | Standardize GFF format
# =================================================================================================


rule main_ref_fix_ids:
    input:
        gff=MAIN_GFF,
        config=rules.agat_config.output,
    output:
        fixed_ID=INT_REFS_DIR / "main_ref_fixed_ID.gff",
    log:
        LOGS / "references" / "annotation" / "main_ref_fix_ids.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gxf2gxf.pl "
        "-g {input.gff} "
        "-o {output.fixed_ID} "
        "-c {input.config} "
        "&> {log}"


rule main_ref_add_locus_tag:
    input:
        fixed_ID=rules.main_ref_fix_ids.output.fixed_ID,
        config=rules.agat_config.output,
    output:
        fixed_locus=INT_REFS_DIR / "main_ref_fixed_locus.gff",
    log:
        LOGS / "references" / "annotation" / "main_ref_add_locus_tag.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sq_add_locus_tag.pl "
        "--gff {input.fixed_ID} "
        "--li ID "
        "-o {output.fixed_locus} "
        "-c {input.config} "
        "&> {log}"


rule main_ref_fix_descriptions:
    input:
        fixed_locus=rules.main_ref_add_locus_tag.output.fixed_locus,
        config=rules.agat_config.output,
    output:
        fixed_description=INT_REFS_DIR / "main_ref_fixed_description.gff",
    log:
        LOGS / "references" / "annotation" / "main_ref_fix_descriptions.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_manage_attributes.pl "
        "--gff {input.fixed_locus} "
        "--tag product/description "
        "-o {output.fixed_description} "
        "-c {input.config} "
        "&> {log}"


rule main_ref_gff2tsv:
    input:
        fixed_description=rules.main_ref_fix_descriptions.output.fixed_description,
        config=rules.agat_config.output,
    output:
        tsv=INT_REFS_DIR / "main_ref_fixed.tsv",
    log:
        LOGS / "references" / "annotation" / "main_ref_gff2tsv.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "--gff {input.fixed_description} "
        "-o {output.tsv} "
        "-c {input.config} "
        "&> {log}"


rule main_ref_recreate_ids:
    input:
        tsv=rules.main_ref_gff2tsv.output.tsv,
    output:
        gff=INT_REFS_DIR / "main_ref.gff",
        tsv=INT_REFS_DIR / "main_ref.tsv",
    log:
        LOGS / "references" / "annotation" / "main_ref_recreate_ids.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "python workflow/scripts/recreate_ids.py "
        "-i {input.tsv} "
        "-og {output.gff} "
        "-ot {output.tsv} "
        "&> {log}"


rule main_ref_symlinks:
    input:
        fasta=MAIN_FASTA,
        gff=rules.main_ref_recreate_ids.output.gff,
    output:
        fasta=INT_REFS_DIR / "{lineage}" / "main_ref.fasta",
        gff=INT_REFS_DIR / "{lineage}" / "main_ref.gff",
    log:
        LOGS / "references" / "annotation" / "main_liks_{lineage}.log",
    conda:
        "../envs/shell.yaml"
    shell:
        "ln -s -r {input.fasta} {output.fasta} &>> {log}"
        "&& "
        "ln -s -r {input.gff} {output.gff} &>> {log}"


# ==================================================================================================
#   Per lineage | Annotate reference genomes using main reference
# ==================================================================================================


rule ref2ref_liftoff:
    input:
        target_refs=INT_REFS_DIR / "{lineage}" / "{lineage}.fasta",
        fasta=rules.main_ref_symlinks.output.fasta,
        gff=rules.main_ref_symlinks.output.gff,
    output:
        target_gff=INT_REFS_DIR / "{lineage}" / "liftoff.gff_polished",
        unmapped=INT_REFS_DIR / "{lineage}" / "unmapped_features.txt",
        intermediate=directory(INT_REFS_DIR / "{lineage}" / "intermediate_liftoff"),
        fai_main=INT_REFS_DIR / "{lineage}" / "main_ref.fasta.fai",
        fai=INT_REFS_DIR / "{lineage}" / "{lineage}.fasta.fai",
        mmi=INT_REFS_DIR / "{lineage}" / "{lineage}.fasta.mmi",
    params:
        refdir=INT_REFS_DIR,
        extra=config["annotate_references"]["liftoff"]["extra"],
    log:
        LOGS / "references" / "annotation" / "ref2ref_liftoff_{lineage}.log",
    threads: config["annotate_references"]["liftoff"]["threads"]
    resources:
        tmpdir=TEMPDIR,
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
        "&> {log}"

rule rename_polished:
    input:
        gff=rules.ref2ref_liftoff.output.target_gff,
    output:
        gff=INT_REFS_DIR / "{lineage}" / "{lineage}_annotated.gff",
    log:
        LOGS / "references" / "annotation" / "rename_polished_{lineage}.log",
    conda:
        "../envs/shell.yaml"
    shell:
        "mv {input.gff} {output.gff} "
        "&> {log}"

rule refs_unmapped_features:
    input:
        DATASET_DIR / "metadata.csv",
        rules.main_ref_recreate_ids.output.tsv,
        expand(rules.ref2ref_liftoff.output.unmapped, lineage=LINEAGES),
    output:
        REFS_DIR / "refs_unmapped_features.tsv",
    conda:
        "../envs/pandas.yaml"
    params:
        refdir=INT_REFS_DIR,
    log:
        LOGS / "references" / "annotation" / "refs_unmapped_features.log",
    script:
        "../scripts/refs_unmapped_features.py"



