# =================================================================================================
#   Join lineages | Create a single table with the annotation of all lineages
# =================================================================================================


rule join_ref_annotations:
    input:
        expand(REFS_DIR / "{lineage}" / "{lineage}.gff.tsv", lineage=LINEAGES),
    output:
        REFS_DIR / "all_lineages.gff.tsv",
    log:
        "logs/references/join_ref_annotations.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_ref_annotations.py"


# =================================================================================================
#   Join dataset | Join tables with sequences, mapq_depth, CNV and variant annotation
# =================================================================================================


rule join_sequences:
    input:
        cds=expand(INT_SAMPLES_DIR / "annotation" / "{sample}" / "cds.csv", sample=SAMPLES),
        prots=expand(INT_SAMPLES_DIR / "annotation" / "{sample}" / "proteins.csv", sample=SAMPLES),
    output:
        sequences=INT_DATASET_DIR / "sequences.csv",
    log:
        "logs/dataset/join_sequences.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_sequences.py"


rule join_mapq_depth:
    input:
        expand(
            SAMPLES_DIR / "depth_quality" / "{sample}" / "mapq_depth_by_feature.tsv", sample=SAMPLES
        ),
    output:
        DATASET_DIR / "depth_quality" / "mapq_depth_by_feature.tsv",
    log:
        "logs/dataset/depth_quality/join_mapq_depth.log",
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_tables.py"


rule join_cnv:
    input:
        expand(SAMPLES_DIR / "cnv" / "{sample}" / "cnv_calls.tsv", sample=SAMPLES),
    output:
        DATASET_DIR / "cnv" / "cnv_calls.tsv",
    log:
        "logs/dataset/cnv/join_cnv_calls.log",
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_tables.py"


rule join_variant_annotation:
    input:
        effects=expand(INT_DATASET_DIR / "snpeff" / "{lineage}_effects.tsv", lineage=LINEAGES),
        variants=expand(INT_DATASET_DIR / "snpeff" / "{lineage}_variants.tsv", lineage=LINEAGES),
        lofs=expand(INT_DATASET_DIR / "snpeff" / "{lineage}_lofs.tsv", lineage=LINEAGES),
        nmds=expand(INT_DATASET_DIR / "snpeff" / "{lineage}_nmds.tsv", lineage=LINEAGES),
        presence=expand(INT_DATASET_DIR / "snpeff" / "{lineage}_presence.tsv", lineage=LINEAGES),
    output:
        effects=DATASET_DIR / "snpeff" / "effects.tsv",
        variants=DATASET_DIR / "snpeff" / "variants.tsv",
        lofs=DATASET_DIR / "snpeff" / "lofs.tsv",
        nmds=DATASET_DIR / "snpeff" / "nmds.tsv",
        presence=DATASET_DIR / "snpeff" / "presence.tsv",
    log:
        "logs/dataset/snpeff/join_variant_annotation.log",
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_variant_annotation.py"

# =================================================================================================
#   Dataset | Create final database
# =================================================================================================


rule complete_db:
    input:
        metadata=rules.quality_filter.output.metadata,
        chrom_names=rules.quality_filter.output.chromosomes,
        cnv=rules.join_cnv.output,
        md=rules.join_mapq_depth.output,
        gffs=rules.join_ref_annotations.output,
        effects=rules.join_variant_annotation.output.effects,
        variants=rules.join_variant_annotation.output.variants,
        presence=rules.join_variant_annotation.output.presence,
        lofs=rules.join_variant_annotation.output.lofs,
        nmds=rules.join_variant_annotation.output.nmds,
        seqs=rules.join_sequences.output,
    output:
        DATASET_DIR / "database.db",
    log:
        "logs/dataset/complete_db.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "xonsh workflow/scripts/build_database.xsh "
        "-m {input.metadata} "
        "-ch {input.chrom_names} "
        "-cnv {input.cnv} "
        "-md {input.md} "
        "-g {input.gffs} "
        "-e {input.effects} "
        "-v {input.variants} "
        "-p {input.presence} "
        "-l {input.lofs} "
        "-n {input.nmds} "
        "-s {input.seqs} "
        "-o {output} &> {log}"
