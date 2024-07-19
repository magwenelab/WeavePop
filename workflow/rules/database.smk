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
#   Join dataset | Add sample sequences to database. (Runs per sample but output is dataset-wide)
# =================================================================================================
# Make SQL database with cds of all samples
rule cds2db:
    input: 
        cds = rules.agat_cds.output.cds,
    output:
        touch(DATASET_OUTDIR / "database" / "{sample}" / "cds.done"),
    params:
        db = DATASET_OUTDIR / "sequences.db"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/sequences/cds2db_{sample}.log"
    shell:
        "python workflow/scripts/build_sequences_db.py "
        "-d {params.db} "
        "-f {input.cds} "
        "-s {wildcards.sample} "
        "-t DNA "
        "&> {log}"

# Make SQL database with proteins of all samples
rule prots2db:
    input: 
        prots = rules.agat_prots.output.prots,
    output:
        touch(DATASET_OUTDIR / "database" / "{sample}" / "prots.done")
    params:
        db = DATASET_OUTDIR / "sequences.db"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/sequences/prots2db_{sample}.log"
    shell:
        "python workflow/scripts/build_sequences_db.py "
        "-d {params.db} "
        "-f {input.prots} "
        "-s {wildcards.sample} "
        "-t PROTEIN "
        "&> {log}"

# =================================================================================================
#   Per dataset | Join feature MAPQ and Depth
# =================================================================================================

rule join_mapq_depth:
    input:
        expand(rules.mapq_depth.output.tsv,sample=SAMPLES)
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
        expand(rules.cnv_calling.output,sample=SAMPLES),
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
rule join_dataframes:
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
        effects = rules.join_dataframes.output.effects,
        variants = rules.join_dataframes.output.variants,
        presence = rules.join_dataframes.output.presence,
        lofs = rules.join_dataframes.output.lofs,
        nmds = rules.join_dataframes.output.nmds,
        cds = expand(DATASET_OUTDIR / "database" / "{sample}" / "cds.done", sample=SAMPLES),
        prots = expand(DATASET_OUTDIR / "database" / "{sample}" / "prots.done", sample=SAMPLES),
    output:
        DATASET_OUTDIR / "database.db"
    params:
        sequences = DATASET_OUTDIR / "sequences.db",
        column = 'group'
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
        "-s {params.sequences} "
        "-c {params.column} "
        "-o {output} &> {log}"