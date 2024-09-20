# =================================================================================================
#   Join lineages | Create a single GFF file with all lineages
# =================================================================================================
rule join_gffs:
    input:
        expand(INT_REFS_DIR / "{lineage}" / "{lineage}.gff.tsv", lineage=LINEAGES)
    output:
        INT_REFS_DIR / "all_lineages.gff.tsv"
    log:
        "logs/references/join_gffs.log"
    resources:
        tmpdir = TEMPDIR
    shell:
        "python workflow/scripts/join_gffs.py -o {output} {input} &> {log}"

# =================================================================================================
#   Join dataset | Join tables with sequences and convert them to SQL db
# =================================================================================================

rule join_sequences:
    input:
        cds = expand(INT_SAMPLES_DIR / "annotation" / "{sample}" / "cds.csv", sample=SAMPLES),
        prots = expand(INT_SAMPLES_DIR / "annotation" / "{sample}" / "proteins.csv", sample=SAMPLES)
    output:
        joined = INT_DATASET_DIR / "sequences.csv"
    log:
        "logs/dataset/join_sequences.log"
    run: 
        with open(log[0], "w") as log_file:
            try:
                cds = pd.concat([pd.read_csv(f) for f in input.cds])
                prots = pd.concat([pd.read_csv(f) for f in input.prots])
                sequences = pd.concat([cds, prots])
                sequences.to_csv(output.joined, sep=",", index=False, header = True)
                log_file.write("Successfully concatenated sequence tables.\n")
            except Exception as e:
                log_file.write(f"Error: {e}\n")
                raise e
# =================================================================================================
#   Per dataset | Join feature MAPQ and Depth
# =================================================================================================

rule join_mapq_depth:
    input:
        expand(SAMPLES_DIR / "depth_quality" / "{sample}" / "feature_mapq_depth.tsv",sample=SAMPLES)
    output:
        DATASET_DIR / "depth_quality" / "feature_mapq_depth.tsv"
    log:
        "logs/dataset/depth_quality/join_mapq_depth.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp
        """

# =================================================================================================
#   Per dataset | Join CNV calls
# =================================================================================================

rule join_cnv_calling:
    input:
        expand(SAMPLES_DIR / "cnv" / "{sample}" / "cnv_calls.tsv",sample=SAMPLES),
    output:
        DATASET_DIR / "cnv" / "cnv_calls.tsv",
    log:
        "logs/dataset/cnv/join_cnv_calls.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp
        """

# =================================================================================================
#   Per dataset | Create final database
# =================================================================================================
# Join all effects, variants, lofs, nmds and presence tables
rule join_variant_annotation:
    input:
        effects = expand(INT_DATASET_DIR / "snps" / "{lineage}_effects.tsv", lineage=LINEAGES),
        variants = expand(INT_DATASET_DIR / "snps" / "{lineage}_variants.tsv", lineage=LINEAGES),
        lofs = expand(INT_DATASET_DIR / "snps" / "{lineage}_lofs.tsv", lineage=LINEAGES),
        nmds = expand(INT_DATASET_DIR / "snps" / "{lineage}_nmds.tsv", lineage=LINEAGES),
        presence = expand(INT_DATASET_DIR / "snps" / "{lineage}_presence.tsv", lineage=LINEAGES)
    output:
        effects = DATASET_DIR / "snps" / "effects.tsv",
        variants = DATASET_DIR / "snps" / "variants.tsv",
        lofs = DATASET_DIR / "snps" / "lofs.tsv",
        nmds = DATASET_DIR / "snps" / "nmds.tsv",
        presence = DATASET_DIR / "snps" / "presence.tsv"
    log:
        "logs/dataset/snps/join_variant_annotation.log"
    run:
        with open(log[0], "w") as log_file:
            try:
                effects = pd.concat([pd.read_csv(f, sep="\t") for f in input.effects])
                variants = pd.concat([pd.read_csv(f, sep="\t") for f in input.variants])
                lofs = pd.concat([pd.read_csv(f, sep="\t") for f in input.lofs])
                nmds = pd.concat([pd.read_csv(f, sep="\t") for f in input.nmds])
                presence = pd.concat([pd.read_csv(f, sep="\t") for f in input.presence])
                effects.to_csv(output.effects, sep = "\t", index=False)
                variants.to_csv(output.variants, sep="\t", index=False)
                lofs.to_csv(output.lofs, sep="\t", index=False)
                nmds.to_csv(output.nmds, sep="\t", index=False)
                presence.to_csv(output.presence, sep="\t", index=False)
                log_file.write("Successfully concatenated variant annotation tables.\n")
            except Exception as e:
                log_file.write(f"Error: {e}\n")
                raise e

# Create the final database
rule complete_db:
    input:
        metadata = INT_DATASET_DIR / "metadata.csv",
        chrom_names = rules.copy_config.output.c,
        cnv = rules.join_cnv_calling.output,
        md = rules.join_mapq_depth.output,
        gffs = rules.join_gffs.output,
        effects = rules.join_variant_annotation.output.effects,
        variants = rules.join_variant_annotation.output.variants,
        presence = rules.join_variant_annotation.output.presence,
        lofs = rules.join_variant_annotation.output.lofs,
        nmds = rules.join_variant_annotation.output.nmds,
        seqs = rules.join_sequences.output
    output:
        DATASET_DIR / "database.db"
    log:
        "logs/dataset/complete_db.log"
    resources:
        tmpdir = TEMPDIR
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