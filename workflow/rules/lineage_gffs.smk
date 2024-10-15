# =================================================================================================
#   Per lineage | Standardize GFF format and convert to TSV
# =================================================================================================


rule agat_convert_gxf2gxf:
    input:
        gff=REF_DATA / "{lineage}.gff",
        config=rules.agat_config.output,
    output:
        fixed_ID=INT_REFS_DIR / "{lineage}" / "{lineage}.fixed.gff",
    log:
        "logs/references/agat_convert_gxf2gxf_{lineage}.log",
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


rule agat_add_locus_tag:
    input:
        fixed_ID=rules.agat_convert_gxf2gxf.output.fixed_ID,
        config=rules.agat_config.output,
    output:
        fixed_locus=INT_REFS_DIR / "{lineage}" / "{lineage}.fixed_locus.gff",
    log:
        "logs/references/agat_add_locus_tag_{lineage}.log",
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


rule agat_manage_attributes:
    input:
        fixed_locus=rules.agat_add_locus_tag.output.fixed_locus,
        config=rules.agat_config.output,
    output:
        fixed_description=INT_REFS_DIR / "{lineage}" / "{lineage}.fixed_description.gff",
    log:
        "logs/references/agat_manage_attributes_{lineage}.log",
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


rule agat_convert_gff2tsv:
    input:
        fixed_description=rules.agat_manage_attributes.output.fixed_description,
        config=rules.agat_config.output,
    output:
        tsv=INT_REFS_DIR / "{lineage}" / "{lineage}.tsv",
    log:
        "logs/references/agat_convert_gff2tsv_{lineage}.log",
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


# Recreate IDs
rule fix_gff_tsv:
    input:
        tsv=rules.agat_convert_gff2tsv.output.tsv,
    output:
        gff=REFS_DIR / "{lineage}.gff",
        tsv=INT_REFS_DIR / "{lineage}" / "{lineage}.gff.tsv",
    log:
        "logs/references/fix_gff_tsv_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    shell:
        "python workflow/scripts/fix_gff.py "
        "-i {input.tsv} "
        "-og {output.gff} "
        "-ot {output.tsv} "
        "&> {log}"

rule gff2bed:
    input:
        rules.agat_convert_gff2tsv.output,
    output:
        INT_REFS_DIR / "{lineage}" / "{lineage}.gff.bed",
    log:
        "logs/references/gff2bed_{lineage}.log",
    conda:
        "../envs/shell.yaml"
    shell:
        """
        tail -n +2 {input} | awk -F '\t' -v OFS='\t' '{{print $1"\t"($4-1)"\t"$5"\t"$13}}'  1> {output} 2> {log}
        """


# Intersect repetitive sequences with genetic features

rule intersect_features_repeats:
    input:
        features=rules.gff2bed.output,
        repeats=REFS_DIR / "{lineage}_repeats.bed",
    output:
        INT_REFS_DIR / "{lineage}" / "feature_repeats.tsv",
    log:
        "logs/references/intersect_features_repeats_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/intersect_features_repeats.xsh "
        "-f {input.features} "
        "-r {input.repeats} "
        "-o {output} "
        "&> {log}"