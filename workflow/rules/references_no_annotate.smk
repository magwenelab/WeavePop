# =================================================================================================
#   Per lineage | Standardize GFF format and convert to TSV
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
        gff=INT_REFS_DIR / "{lineage}" / "{lineage}_2.gff",
        tsv=INT_REFS_DIR / "{lineage}" / "{lineage}_2.gff.tsv",
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


# ==================================================================================================
#   Per lineage | Add repetitive sequences, introns, intergenic regions, and convert to TSV
# ==================================================================================================


rule ref_add_intergenic:
    input:
        gff=rules.ref_recreate_ids.output.gff,
        config=rules.agat_config.output,
    output:
        INT_REFS_DIR / "{lineage}" / "{lineage}_intergenic.gff",
    params:
        extra=config["agat"]["extra"],
    log:
        LOGS / "references" / "no_annotation" / "ref_add_intergenic_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_add_intergenic_regions.pl "
        "-g {input.gff} "
        "-o {output} "
        "-c {input.config} "
        "{params.extra} "
        "&> {log}"

rule ref_add_introns:
    input:
        gff=rules.ref_add_intergenic.output,
        config=rules.agat_config.output,
    output:
        gff=INT_REFS_DIR / "{lineage}" / "{lineage}_interg_introns.gff",
    params:
        extra=config["agat"]["extra"],
    log:
        LOGS / "references" / "no_annotation" / "ref_add_introns_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_add_introns.pl "
        "-g {input.gff} "
        "-o {output.gff} "
        "-c {input.config} "
        "{params.extra} "
        "&> {log}"


rule ref_gff2tsv:
    input:
        target=rules.ref_add_introns.output.gff,
        config=rules.agat_config.output,
    output:
        tsv=INT_REFS_DIR / "{lineage}" / "{lineage}_interg_introns.gff.tsv",
    log:
        LOGS / "references" / "no_annotation" / "gff2tsv_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "-gff {input.target} "
        "-c {input.config} "
        "-o {output.tsv} "
        "&> {log}"

rule ref_add_repeats:
    input:
        gff=rules.ref_add_introns.output.gff,
        repeats=REFS_DIR / "{lineage}" / "{lineage}_repeats.bed",
    output:
        INT_REFS_DIR / "{lineage}" / "{lineage}_repeats.gff",
    log:
        LOGS / "references" / "no_annotation" / "ref_add_repeats_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/ref_add_repeats.xsh "
        "-g {input.gff} "
        "-r {input.repeats} "
        "-o {output} "
        "&> {log}"


rule ref_gff2tsv_2:
    input:
        target=rules.ref_add_repeats.output,
        config=rules.agat_config.output,
    output:
        tsv=INT_REFS_DIR / "{lineage}" / "{lineage}_repeats.gff.tsv",
    log:
        LOGS / "references" / "no_annotation" / "gff2tsv_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "-gff {input.target} "
        "-c {input.config} "
        "-o {output.tsv} "
        "&> {log}"

rule ref_reformat_annotation:
    input:
        tsv=rules.ref_gff2tsv_2.output.tsv,
    output:
        tsv=REFS_DIR / "{lineage}" / "{lineage}.gff.tsv",
        gff=REFS_DIR / "{lineage}" / "{lineage}.gff",
    log:
        LOGS / "references" / "no_annotation" / "ref_reformat_annotation_{lineage}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/reformat_annotation.py"





