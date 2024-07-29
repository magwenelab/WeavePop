# =================================================================================================
#   Join lineages | Create a single GFF file with all lineages
# =================================================================================================
rule join_gffs:
    input:
        expand(REFDIR / "{lineage}" / "{lineage}.gff.tsv", lineage=LINEAGES)
    output:
        REFDIR / "all.gff.tsv"
    log:
        "logs/references/join_gffs.log"
    shell:
        "python workflow/scripts/join_gffs.py -o {output} {input} &> {log}"

# =================================================================================================
#   Join dataset | Join tables with sequences and convert them to SQL db
# =================================================================================================

rule join_sequences:
    input:
        cds = expand(OUTDIR / "agat" / "{sample}" / "cds.csv", sample=SAMPLES),
        prots = expand(OUTDIR / "agat" / "{sample}" / "proteins.csv", sample=SAMPLES)
    output:
        joined = DATASET_OUTDIR / "sequences.csv"
    run: 
        cds = pd.concat([pd.read_csv(f) for f in input.cds])
        prots = pd.concat([pd.read_csv(f) for f in input.prots])
        sequences = pd.concat([cds, prots])
        sequences.to_csv(output.joined, sep=",", index=False, header = True)
    
# =================================================================================================
#   Per dataset | Join feature MAPQ and Depth
# =================================================================================================

rule join_mapq_depth:
    input:
        expand(OUTDIR / "samtools" / "{sample}" / "feature_mapq_depth.tsv",sample=SAMPLES)
    output:
        DATASET_OUTDIR / "files" / "feature_mapq_depth.tsv"
    log:
        "logs/dataset/files/join_mapq_depth.log"
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
        expand(OUTDIR / "cnv" / "{sample}" / "cnv_calls.tsv",sample=SAMPLES),
    output:
        DATASET_OUTDIR / "files" / "cnv_calls.tsv",
    log:
        "logs/dataset/files/join_cnv_calls.log"
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
        effects = expand(DATASET_OUTDIR / "snps" / "{lineage}_effects.tsv", lineage=LINEAGES),
        variants = expand(DATASET_OUTDIR / "snps" / "{lineage}_variants.tsv", lineage=LINEAGES),
        lofs = expand(DATASET_OUTDIR / "snps" / "{lineage}_lofs.tsv", lineage=LINEAGES),
        nmds = expand(DATASET_OUTDIR / "snps" / "{lineage}_nmds.tsv", lineage=LINEAGES),
        presence = expand(DATASET_OUTDIR / "snps" / "{lineage}_presence.tsv", lineage=LINEAGES)
    output:
        effects = DATASET_OUTDIR / "snps" / "all_effects.tsv",
        variants = DATASET_OUTDIR / "snps" / "all_variants.tsv",
        lofs = DATASET_OUTDIR / "snps" / "all_lofs.tsv",
        nmds = DATASET_OUTDIR / "snps" / "all_nmds.tsv",
        presence = DATASET_OUTDIR / "snps" / "all_presence.tsv"
    run:
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

# Create the final database
rule complete_db:
    input:
        metadata = SAMPLEFILE,
        chrom_names = CHROM_NAMES,
        sv = rules.join_cnv_calling.output,
        md = rules.join_mapq_depth.output,
        gffs = rules.join_gffs.output,
        effects = rules.join_variant_annotation.output.effects,
        variants = rules.join_variant_annotation.output.variants,
        presence = rules.join_variant_annotation.output.presence,
        lofs = rules.join_variant_annotation.output.lofs,
        nmds = rules.join_variant_annotation.output.nmds,
        seqs = rules.join_sequences.output
    output:
        DATASET_OUTDIR / "database.db"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/complete_db.log"
    shell:
        "xonsh workflow/scripts/build_database.xsh "
        "-m {input.metadata} "
        "-ch {input.chrom_names} "
        "-sv {input.sv} "
        "-md {input.md} "
        "-g {input.gffs} "
        "-e {input.effects} "
        "-v {input.variants} "
        "-p {input.presence} "
        "-l {input.lofs} "
        "-n {input.nmds} "
        "-s {input.seqs} "
        "-o {output} &> {log}"