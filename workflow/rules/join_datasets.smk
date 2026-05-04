# =================================================================================================
#   Join the metadata, chromosomes and mapping stats files from the different datasets
# =================================================================================================


rule join_metadata:
    input:
        expand(os.path.join("{dir}", DATASET_DIR_NAME, "metadata.csv"), dir=LIST_PATHS),
    output:
        DATASET_DIR / "metadata.csv",
    params:
        config["datasets"]["names"].split(","),
    log:
        LOGS / "join_datasets" / "join_metadata.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_metadata.py"


rule join_chromosomes:
    input:
        expand(
            os.path.join("{dir}", DATASET_DIR_NAME, "chromosomes.csv"), dir=LIST_PATHS
        ),
    output:
        DATASET_DIR / "chromosomes.csv",
    log:
        LOGS / "join_datasets" / "join_chromosomes.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_chromosomes.py"


rule join_mapping_stats:
    input:
        expand(
            os.path.join("{dir}", DATASET_DIR_NAME, "depth_quality" , "mapping_stats.tsv"), dir=LIST_PATHS
        ),
    output:
        DATASET_DIR / "depth_quality" / "mapping_stats.tsv",
    log:
        LOGS / "join_datasets" / "join_mapping_stats.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_tables.py"


# =================================================================================================
#   Get names of ref_genomes in combined datasets
# =================================================================================================


checkpoint get_ref_genomes:
    input:
        rules.join_metadata.output,
    output:
        directory(INT_REFS_DIR / "ref_genome_names"),
    log:
        LOGS / "join_datasets" / "get_ref_genomes.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/get_ref_genomes.py"


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
        "../envs/pandas.yaml"
    log:
        LOGS / "join_datasets" / "join_sequences.log",
    script:
        "../scripts/join_sequences.py"


rule join_cnv:
    input:
        input_join_cnv,
    output:
        DATASET_DIR / "cnv" / "cnv_calls.tsv",
    log:
        LOGS / "join_datasets" / "join_cnv.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_tables.py"


rule join_cnv_chromosomes:
    input:
        input_join_cnv_chrom,
    output:
        DATASET_DIR / "cnv" / "cnv_chromosomes.tsv",
    log:
        LOGS / "join_datasets" / "join_cnv_chromosomes.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_tables.py"


rule join_mapq_depth:
    input:
        input_join_mapq_depth,
    output:
        DATASET_DIR / "depth_quality" / "mapq_depth_by_feature.tsv",
    log:
        LOGS / "join_datasets" / "join_mapq_depth.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_tables.py"


rule join_ref_annotations:
    input:
        input_join_ref_annotations,
    output:
        INT_REFS_DIR / "all_ref_genomes_gff.tsv",
    log:
        LOGS / "join_datasets" / "join_ref_annotations.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_ref_annotations.py"


rule join_ref_sequences:
    input:
        cds=input_join_ref_cds,
        prots=input_join_ref_prots,
    output:
        sequences=INT_REFS_DIR / "all_refs_sequences.csv",
    log:
        LOGS / "join_datasets" / "join_ref_sequences.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_ref_sequences.py"


# =================================================================================================
#   Redo union of VCFs and variant annotation
# =================================================================================================


rule copy_snpeff_data:
    input:
        input_copy_speff_data,
    output:
        INT_REFS_DIR / "snpeff_data" / "copy.done",
    params:
        dir=INT_REFS_DIR / "snpeff_data",
    log:
        LOGS / "join_datasets" / "copy_snpeff_data.log",
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
            os.path.join(
                "{dir}", INTDIR_NAME, REFS_DIR_NAME, "snpeff_data", "snpEff.config"
            ),
            dir=LIST_PATHS,
        ),
    output:
        INT_REFS_DIR / "snpeff_data" / "snpEff.config",
    log:
        LOGS / "join_datasets" / "copy_snpeff_config.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/shell.yaml"
    shell:
        "cat {input.config} 1> {output} 2> {log}"


rule unite_vcfs:
    input:
        unpack(input_unite_vcfs),
    output:
        vcf=INT_DATASET_DIR / "snpeff" / "{ref_genome}_union.vcf",
        tsv=INT_DATASET_DIR / "snpeff" / "{ref_genome}_presence.tsv",
    params:
        tmp_dir=os.path.join(TEMPDIR, "tmp_{ref_genome}"),
    log:
        LOGS / "join_datasets" / "unite_vcfs_{ref_genome}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    script:
        "../scripts/unite_vcfs.xsh"


rule snpeff:
    input:
        vcf=rules.unite_vcfs.output.vcf,
        config=rules.copy_snpeff_config.output,
    output:
        vcf=INT_DATASET_DIR / "snpeff" / "{ref_genome}_snpeff.vcf",
        html=INT_DATASET_DIR / "snpeff" / "{ref_genome}_snpeff.html",
    params:
        dir=os.getcwd() / INT_REFS_DIR / "snpeff_data",
        name="Species_name_{ref_genome}",
    log:
        LOGS / "join_datasets" / "snpeff_{ref_genome}.log",
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
        tsv=INT_REFS_DIR / "{ref_genome}.gff.tsv",
    log:
        LOGS / "join_datasets" / "symlink_ref_gff_{ref_genome}.log",
    conda:
        "../envs/shell.yaml"
    shell:
        "ln -sr {input} {output} &> {log}"


rule extract_vcf_annotation:
    input:
        vcf=rules.snpeff.output.vcf,
        tsv=rules.symlink_ref_gff.output.tsv,
        metadata=DATASET_DIR / "metadata.csv",
        presence=rules.unite_vcfs.output.tsv,
    output:
        effects=INT_DATASET_DIR / "snpeff" / "{ref_genome}_effects.tsv",
        variants=INT_DATASET_DIR / "snpeff" / "{ref_genome}_variants.tsv",
        lofs=INT_DATASET_DIR / "snpeff" / "{ref_genome}_lofs.tsv",
        nmds=INT_DATASET_DIR / "snpeff" / "{ref_genome}_nmds.tsv",
        classif=INT_DATASET_DIR / "snpeff" / "{ref_genome}_variant_classification.tsv",
    log:
        LOGS / "join_datasets" / "extract_vcf_annotation_{ref_genome}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    script:
        "../scripts/extract_vcf_annotation.py"


rule join_variant_annotation:
    input:
        effects=expand(
            INT_DATASET_DIR / "snpeff" / "{ref_genome}_effects.tsv",
            ref_genome=REF_GENOMES,
        ),
        variants=expand(
            INT_DATASET_DIR / "snpeff" / "{ref_genome}_variants.tsv",
            ref_genome=REF_GENOMES,
        ),
        classif=expand(
            INT_DATASET_DIR / "snpeff" / "{ref_genome}_variant_classification.tsv",
            ref_genome=REF_GENOMES,
        ),
        lofs=expand(
            INT_DATASET_DIR / "snpeff" / "{ref_genome}_lofs.tsv",
            ref_genome=REF_GENOMES,
        ),
        nmds=expand(
            INT_DATASET_DIR / "snpeff" / "{ref_genome}_nmds.tsv",
            ref_genome=REF_GENOMES,
        ),
        presence=expand(
            INT_DATASET_DIR / "snpeff" / "{ref_genome}_presence.tsv",
            ref_genome=REF_GENOMES,
        ),
    output:
        effects=DATASET_DIR / "snpeff" / "effects.tsv",
        variants=DATASET_DIR / "snpeff" / "variants.tsv",
        classif=DATASET_DIR / "snpeff" / "variant_classification.tsv",
        lofs=DATASET_DIR / "snpeff" / "lofs.tsv",
        nmds=DATASET_DIR / "snpeff" / "nmds.tsv",
        presence=DATASET_DIR / "snpeff" / "presence.tsv",
    log:
        LOGS / "join_datasets" / "join_variant_annotation.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_variant_annotation.py"


# =================================================================================================
#   Create final database
# =================================================================================================


rule database:
    input:
        metadata=DATASET_DIR / "metadata.csv",
        chrom_names=DATASET_DIR / "chromosomes.csv",
        cnv=rules.join_cnv.output,
        cnv_chromosomes=rules.join_cnv_chromosomes.output,
        md=rules.join_mapq_depth.output,
        gffs=rules.join_ref_annotations.output,
        effects=rules.join_variant_annotation.output.effects,
        variants=rules.join_variant_annotation.output.variants,
        classif=rules.join_variant_annotation.output.classif,
        presence=rules.join_variant_annotation.output.presence,
        lofs=rules.join_variant_annotation.output.lofs,
        nmds=rules.join_variant_annotation.output.nmds,
        seqs=rules.join_sequences.output.sequences,
        ref_seqs=rules.join_ref_sequences.output.sequences,
    output:
        db=DATASET_DIR / "database.db",
    log:
        LOGS / "join_datasets" / "complete_db.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    script:
        "../scripts/database.py"


# =================================================================================================
#   Create plots
# =================================================================================================

rule dataset_summary_plot:
    input:
        metadata=rules.join_metadata.output,
        chroms=rules.join_chromosomes.output,
        stats=rules.join_mapping_stats.output,
    output:
        plot=DATASET_DIR / "plots" / "mapping_summary.png",
    log:
        LOGS / "dataset" / "plots" / "dataset_summary.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_summary_plot.R"

rule dataset_depth_plot:
    input:
        metadata=rules.join_metadata.output,
        chroms=rules.join_chromosomes.output,
        cnv=rules.join_cnv_chromosomes.output,
    output:
        plot=DATASET_DIR / "plots" / "depth_summary.png",
    params:
        column=COLOR_BY,
    log:
        LOGS / "dataset" / "plots" / "dataset_depth_plot.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_depth_plot.R"


rule dataset_variants_plot:
    input:
        metadata=rules.join_metadata.output,
        classif=DATASET_DIR / "snpeff" / "variant_classification.tsv",
        variants=DATASET_DIR / "snpeff" / "variants.tsv",
        presence=DATASET_DIR / "snpeff" / "presence.tsv",
    output:
        plot=DATASET_DIR / "plots" / "variant_summary.png",
    log:
        LOGS
        / "dataset"
        / "plots"
        / "dataset_variants_plot.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_variants_plot.R"

rule dataset_cnv_plot:
    input:
        metadata=rules.join_metadata.output,
        chromosomes=rules.join_chromosomes.output,
        cnv=DATASET_DIR / "cnv" / "cnv_chromosomes.tsv",
    output:
        plot=DATASET_DIR / "plots" / "cnv_summary.png",
    log:
        LOGS
        / "dataset"
        / "plots"
        / "dataset_cnv_plot.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_cnv_plot.R"
