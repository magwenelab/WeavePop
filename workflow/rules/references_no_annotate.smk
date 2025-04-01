# =================================================================================================
#   Per lineage | Standardize GFF format
# =================================================================================================


rule ref_fix_ids:
    input:
        gff=REF_DATA / "{lineage}.gff",
        config=rules.agat_config.output,
    output:
        fixed_ID=INT_REFS_DIR / "{lineage}" / "{lineage}.fixed.gff",
    log:
        LOGS / "references" / "no_annotation" / "ref_fix_ids_{lineage}.log",
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


rule ref_add_locus_tag:
    input:
        fixed_ID=rules.ref_fix_ids.output.fixed_ID,
        config=rules.agat_config.output,
    output:
        fixed_locus=INT_REFS_DIR / "{lineage}" / "{lineage}.fixed_locus.gff",
    log:
        LOGS / "references" / "no_annotation" / "ref_add_locus_tag_{lineage}.log",
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


rule ref_fix_descriptions:
    input:
        fixed_locus=rules.ref_add_locus_tag.output.fixed_locus,
        config=rules.agat_config.output,
    output:
        fixed_description=INT_REFS_DIR / "{lineage}" / "{lineage}.fixed_description.gff",
    log:
        LOGS / "references" / "no_annotation" / "ref_fix_descriptions_{lineage}.log",
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


rule ref_gff2tsv_1:
    input:
        gff=rules.ref_fix_descriptions.output.fixed_description,
        config=rules.agat_config.output,
    output:
        tsv=INT_REFS_DIR / "{lineage}" / "{lineage}_1.gff.tsv",
    log:
        LOGS / "references" / "no_annotation" / "ref_gff2tsv_1_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "--gff {input.gff} "
        "-o {output.tsv} "
        "-c {input.config} "
        "&> {log}"


rule ref_recreate_ids:
    input:
        tsv=rules.ref_gff2tsv_1.output.tsv,
    output:
        gff=INT_REFS_DIR / "{lineage}" / "{lineage}_annotated.gff",
        tsv=INT_REFS_DIR / "{lineage}" / "{lineage}_annotated.gff.tsv",
    log:
        LOGS / "references" / "no_annotation" / "ref_recreate_ids_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "python workflow/scripts/recreate_ids.py "
        "-i {input.tsv} "
        "-og {output.gff} "
        "-ot {output.tsv} "
        "&> {log}"




