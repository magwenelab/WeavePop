# Join all lineages gffs into one
rule join_gffs:
    input:
        expand(REFDIR / "{lineage}" / "{lineage}.gff.tsv", lineage=LINEAGES)
    output:
        REFDIR / "all.gff.tsv"
    log:
        "logs/references/join_gffs.log"
    shell:
        "python workflow/scripts/join_gffs.py -o {output} {input} &> {log}"

# Extract cds and protein sequences from reference genomes
rule extract_ref_seqs:
    input:
        gff = REFDIR / "{lineage}" / "{lineage}.gff",
        fasta = REFDIR / "{lineage}" / "{lineage}.fasta",
        config = rules.agat_config.output
    output:
        cds = REFDIR / "{lineage}" / "{lineage}.cds.fa",
        prots = REFDIR / "{lineage}" / "{lineage}.prots.fa"
    conda:
        "../envs/agat.yaml"
    log:
        "logs/references/extract_ref_seqs_{lineage}.log"
    shell:
        "agat_sp_extract_sequences.pl -g {input.gff} -f {input.fasta} -o {output.cds} -c {input.config} &> {log} "
        "&& "
        "agat_sp_extract_sequences.pl -g {input.gff} -f {input.fasta} -o {output.prots} -p -c {input.config} &>> {log} "

# Make symbolic links in the snpeff_data directory and create config file
rule prepare_refs_db:
    input:
        gff = rules.extract_ref_seqs.input.gff,
        fasta = rules.extract_ref_seqs.input.fasta,
        cds = rules.extract_ref_seqs.output.cds,
        prots = rules.extract_ref_seqs.output.prots
    output:
        gff = REFDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "genes.gff",
        fasta = REFDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "sequences.fa",
        cds = REFDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "cds.fa",
        prots = REFDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "protein.fa"
    conda:
        "../envs/variants.yaml"
    params:
        name = config["species_name"] + "_{lineage}",
        config = REFDIR / "snpeff_data" / "snpEff.config",
    log:
        "logs/references/prepare_dbs_{lineage}.log"
    shell:
        """
        echo "{params.name}.genome : {params.name}" >> {params.config} 2> {log} && 
        ln -s -r {input.gff} {output.gff} &>> {log} && 
        ln -s -r {input.fasta} {output.fasta} &>> {log} && 
        ln -s -r {input.cds} {output.cds} &>> {log} && 
        ln -s -r {input.prots} {output.prots} &>> {log} 
        """

# Build snpeff database for the reference genomes
rule build_refs_db:
    input:
        gff = rules.prepare_refs_db.output.gff,
        fasta = rules.prepare_refs_db.output.fasta,
        cds = rules.prepare_refs_db.output.cds,
        prots = rules.prepare_refs_db.output.prots
    output:
        touch(REFDIR / "snpeff_data" / "{lineage}.done")
    conda:
        "../envs/variants.yaml"
    params:
        config = REFDIR/ "snpeff_data" / "snpEff.config",
        dir = os.getcwd() / REFDIR / "snpeff_data",
        name = config["species_name"] + "_{lineage}"
    log:
        "logs/references/build_dbs_{lineage}.log"
    shell:
        """
        snpEff build -gff3 -v -dataDir {params.dir} -config {params.config} {params.name} &>> {log}
        """

def intersect_vcfs_input(wildcards):
    l = LINEAGE_REFERENCE.loc[wildcards.lineage,]
    return {
        "vcfs" : expand(rules.snippy.output.vcf, sample=l["sample"])
    }
    
rule intersect_vcfs:
    input:
        unpack(intersect_vcfs_input)
    output:
        vcf = DATASET_OUTDIR / "snps" / "{lineage}_intersection.vcf",
        tsv = DATASET_OUTDIR / "snps" / "{lineage}_presence.tsv"
    params:
        tmp_dir = DATASET_OUTDIR / 'tmp_{lineage}'
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/intersect_vcfs_{lineage}.log"
    shell:
        "xonsh workflow/scripts/intersect_vcfs.xsh -v {output.vcf} -p {output.tsv} -l {wildcards.lineage} -t {params.tmp_dir} {input.vcfs} &> {log}"

rule snpeff:
    input:
        vcf = rules.intersect_vcfs.output.vcf,
        db_done = rules.build_refs_db.output
    output:
        vcf = DATASET_OUTDIR / "snps" / "{lineage}_snpeff.vcf",
        html = DATASET_OUTDIR / "snps" / "{lineage}_snpeff.html"
    params:
       dir = os.getcwd() / REFDIR / "snpeff_data",
       config = REFDIR / "snpeff_data" / "snpEff.config",
       name = config["species_name"] + "_{lineage}"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/snpeff_{lineage}.log"
    shell:
        "snpEff ann -v -classic -dataDir {params.dir} -config {params.config} -s {output.html} {params.name} {input.vcf} 1> {output.vcf} 2> {log}"

rule extract_vcf_annotation:
    input:
        vcf = rules.snpeff.output.vcf
    output:
        effects = temp(DATASET_OUTDIR / "snps" / "{lineage}_effects.tsv"),
        variants = temp(DATASET_OUTDIR / "snps" / "{lineage}_variants.tsv"),
        lofs = temp(DATASET_OUTDIR / "snps" / "{lineage}_lofs.tsv"),
        nmds = temp(DATASET_OUTDIR / "snps" / "{lineage}_nmds.tsv")
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/extract_vcf_annotation_{lineage}.log"
    shell:
        "xonsh workflow/scripts/extract_vcf_annotation.xsh -i {input.vcf} -e {output.effects} -v {output.variants} -f {output.lofs} -n {output.nmds} -l {wildcards.lineage} &> {log}"

rule join_dataframes:
    input:
        effects = expand(DATASET_OUTDIR / "snps" / "{lineage}_effects.tsv", lineage=LINEAGES),
        variants = expand(DATASET_OUTDIR / "snps" / "{lineage}_variants.tsv", lineage=LINEAGES),
        lofs = expand(DATASET_OUTDIR / "snps" / "{lineage}_lofs.tsv", lineage=LINEAGES),
        nmds = expand(DATASET_OUTDIR / "snps" / "{lineage}_nmds.tsv", lineage=LINEAGES),
        presence = expand(DATASET_OUTDIR / "snps" / "{lineage}_presence.tsv", lineage=LINEAGES)
    output:
        effects = temp(DATASET_OUTDIR / "snps" / "all_effects.tsv"),
        variants = temp(DATASET_OUTDIR / "snps" / "all_variants.tsv"),
        lofs = temp(DATASET_OUTDIR / "snps" / "all_lofs.tsv"),
        nmds = temp(DATASET_OUTDIR / "snps" / "all_nmds.tsv"),
        presence = temp(DATASET_OUTDIR / "snps" / "all_presence.tsv")
    run:
        effects = pd.concat([pd.read_csv(f, sep="\t") for f in input.effects])
        variants = pd.concat([pd.read_csv(f, sep="\t") for f in input.variants])
        lofs = pd.concat([pd.read_csv(f, sep="\t") for f in input.lofs])
        nmds = pd.concat([pd.read_csv(f, sep="\t") for f in input.nmds])
        presence = pd.concat([pd.read_csv(f, sep="\t") for f in input.presence])
        effects.to_csv(output.effects, sep="\t", index=False)
        variants.to_csv(output.variants, sep="\t", index=False)
        lofs.to_csv(output.lofs, sep="\t", index=False)
        nmds.to_csv(output.nmds, sep="\t", index=False)
        presence.to_csv(output.presence, sep="\t", index=False)

# Create the final database
rule complete_db:
    input:
        metadata = SAMPLEFILE,
        chrom_names = CHROM_NAMES,
        sv = rules.dataset_metrics.output.allsv,
        mc = rules.dataset_metrics.output.allmc,
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
        "-mc {input.mc} "
        "-g {input.gffs} "
        "-e {input.effects} "
        "-v {input.variants} "
        "-p {input.presence} "
        "-l {input.lofs} "
        "-n {input.nmds} "
        "-s {params.sequences} "
        "-c {params.column} "
        "-o {output} &> {log}"