# =================================================================================================
#   Join the metadata and chromosomes files from the different datasets
# =================================================================================================


rule join_metadata:
    input:
        expand(os.path.join("{dir}", DATASET_DIR_NAME, "metadata.csv"), dir=LIST_PATHS),
    output:
        DATASET_DIR / "metadata.csv",
    params:
        config["datasets_names"].split(","),
    log:
        "logs/join_datasets/join_metadata.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_metadata.py"


rule join_chromosomes:
    input:
        expand(os.path.join("{dir}", INTDIR_NAME, REFS_DIR_NAME, "chromosomes.csv"), dir=LIST_PATHS),
    output:
        INT_REFS_DIR / "chromosomes.csv",
    log:
        "logs/join_datasets/join_chromosomes.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_chromosomes.py"


checkpoint get_lineages:
    input:
        rules.join_metadata.output,
    output:
        directory(INT_REFS_DIR / "lineage_names"),
    log:
        "logs/join_datasets/get_lineages.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/get_lineages.py"


# =================================================================================================
#   Join tabular results
# =================================================================================================


rule join_sequences:
    input:
        cds=input_join_cds,
        prots=input_join_prots,
    output:
        sequences=INT_DATASET_DIR / "sequences.csv",
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/join_datasets/join_sequences.log",
    script:
        "../scripts/join_sequences.py"


rule join_cnv:
    input:
        input_join_cnv,
    output:
        DATASET_DIR / "cnv" / "cnv_calls.tsv",
    log:
        "logs/join_datasets/join_cnv.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_tables.py"


rule join_mapq_depth:
    input:
        input_join_mapq_depth,
    output:
        DATASET_DIR / "depth_quality" / "feature_mapq_depth.tsv",
    log:
        "logs/join_datasets/join_mapq_depth.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_tables.py"


rule join_ref_annotations:
    input:
        input_join_ref_annotations,
    output:
        INT_REFS_DIR / "all_lineages.gff.tsv",
    log:
        "logs/join_datasets/join_ref_annotations.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_ref_annotations.py"


# =================================================================================================
#   Copy SnpEff lineage databases and config file
# =================================================================================================


rule copy_snpeff_data:
    input:
        input_copy_speff_data,
    output:
        INT_REFS_DIR / "snpeff_data" / "copy.done",
    params:
        dir=INT_REFS_DIR / "snpeff_data",
    log:
        "logs/join_datasets/copy_snpeff_data.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/shell.yaml"
    shell:
        """
        ln -srf {input} {params.dir} 2> {log} && 
        touch {output} 2>> {log}
        """


rule copy_snpeff_config:
    input:
        data=rules.copy_snpeff_data.output,
        config=expand(
            os.path.join("{dir}", INTDIR_NAME, REFS_DIR_NAME, "snpeff_data", "snpEff.config"),
            dir=LIST_PATHS,
        ),
    output:
        INT_REFS_DIR / "snpeff_data" / "snpEff.config",
    log:
        "logs/join_datasets/copy_snpeff_config.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/shell.yaml"
    shell:
        "cat {input.config} 1> {output} 2> {log}"


# =================================================================================================
#   Intersect SNPs, run SnpEff and extract annotations
# =================================================================================================


rule intersect_vcfs:
    input:
        unpack(input_intersect_vcfs),
    output:
        vcf=INT_DATASET_DIR / "snps" / "{lineage}_intersection.vcf",
        tsv=INT_DATASET_DIR / "snps" / "{lineage}_presence.tsv",
    params:
        tmp_dir=os.path.join(TEMPDIR, "tmp_{lineage}"),
    log:
        "logs/join_datasets/intersect_vcfs_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "xonsh workflow/scripts/intersect_vcfs.xsh "
        "-v {output.vcf} "
        "-p {output.tsv} "
        "-l {wildcards.lineage} "
        "-t {params.tmp_dir} "
        "{input.vcfs} "
        "&> {log}"


rule snpeff:
    input:
        vcf=rules.intersect_vcfs.output.vcf,
        config=rules.copy_snpeff_config.output,
    output:
        vcf=INT_DATASET_DIR / "snps" / "{lineage}_snpeff.vcf",
        html=INT_DATASET_DIR / "snps" / "{lineage}_snpeff.html",
    params:
        dir=os.getcwd() / INT_REFS_DIR / "snpeff_data",
        name=config["species_name"] + "_{lineage}",
    log:
        "logs/join_datasets/snpeff_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "snpEff ann "
        "-v "
        "-classic "
        "-dataDir {params.dir} "
        "-config {input.config} "
        "-s {output.html} "
        "{params.name} "
        "{input.vcf} "
        "1> {output.vcf} "
        "2> {log}"


rule symlink_ref_gff:
    input:
        input_symlink_ref_gff,
    output:
        INT_REFS_DIR / "{lineage}.gff.tsv",
    log:
        "logs/join_datasets/symlink_ref_gff_{lineage}.log",
    shell:
        "ln -sr {input} {output} &> {log}"

rule extract_vcf_annotation:
    input:
        vcf=rules.snpeff.output.vcf,
        gff=rules.symlink_ref_gff.output,
    output:
        effects=INT_DATASET_DIR / "snps" / "{lineage}_effects.tsv",
        variants=INT_DATASET_DIR / "snps" / "{lineage}_variants.tsv",
        lofs=INT_DATASET_DIR / "snps" / "{lineage}_lofs.tsv",
        nmds=INT_DATASET_DIR / "snps" / "{lineage}_nmds.tsv",
    log:
        "logs/join_datasets/extract_vcf_annotation_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "xonsh workflow/scripts/extract_vcf_annotation.xsh "
        "-i {input.vcf} "
        "-g {input.gff} "
        "-e {output.effects} "
        "-v {output.variants} "
        "-f {output.lofs} "
        "-n {output.nmds} "
        "-l {wildcards.lineage} "
        "&> {log}"


# =================================================================================================
#   Create final database
# =================================================================================================
rule join_variant_annotation:
    input:
        effects=expand(INT_DATASET_DIR / "snps" / "{lineage}_effects.tsv", lineage=LINEAGES),
        variants=expand(INT_DATASET_DIR / "snps" / "{lineage}_variants.tsv", lineage=LINEAGES),
        lofs=expand(INT_DATASET_DIR / "snps" / "{lineage}_lofs.tsv", lineage=LINEAGES),
        nmds=expand(INT_DATASET_DIR / "snps" / "{lineage}_nmds.tsv", lineage=LINEAGES),
        presence=expand(INT_DATASET_DIR / "snps" / "{lineage}_presence.tsv", lineage=LINEAGES),
    output:
        effects=DATASET_DIR / "snps" / "effects.tsv",
        variants=DATASET_DIR / "snps" / "variants.tsv",
        lofs=DATASET_DIR / "snps" / "lofs.tsv",
        nmds=DATASET_DIR / "snps" / "nmds.tsv",
        presence=DATASET_DIR / "snps" / "presence.tsv",
    log:
        "logs/join_datasets/join_variant_annotation.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/snakemake.yaml"
    script:
        "../scripts/join_variant_annotation.py"


rule complete_db:
    input:
        metadata=DATASET_DIR / "metadata.csv",
        chrom_names=INT_REFS_DIR / "chromosomes.csv",
        cnv=rules.join_cnv.output,
        md=rules.join_mapq_depth.output,
        gffs=rules.join_ref_annotations.output,
        effects=rules.join_variant_annotation.output.effects,
        variants=rules.join_variant_annotation.output.variants,
        presence=rules.join_variant_annotation.output.presence,
        lofs=rules.join_variant_annotation.output.lofs,
        nmds=rules.join_variant_annotation.output.nmds,
        seqs=rules.join_sequences.output.sequences,
    output:
        DATASET_DIR / "database.db",
    log:
        "logs/join_datasets/complete_db.log",
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
        "-o {output} "
        "&> {log}"
