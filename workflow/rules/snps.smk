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


rule complete_db:
    input:
        metadata = SAMPLEFILE,
        chrom_names = CHROM_NAMES,
        sv = rules.dataset_metrics.output.allsv,
        mc = rules.dataset_metrics.output.allmc,
        vcfs = expand(rules.snippy.output.vcf, sample=SAMPLES),
        databases = expand(rules.build_refs_db.output, lineage=LINEAGES),
        gffs = rules.join_gffs.output,
        cds = expand(DATASET_OUTDIR / "database" / "{sample}" / "cds.done", sample=SAMPLES),
        prots = expand(DATASET_OUTDIR / "database" / "{sample}" / "prots.done", sample=SAMPLES)
    output:
        DATASET_OUTDIR / "database.db"
    params:
        column = 'group',
        sp = config["species_name"],
        temp_dir = DATASET_OUTDIR / 'tmp',
        dir = os.getcwd() / REFDIR/ "snpeff_data",
        config = REFDIR / "snpeff_data" / "snpEff.config",
        sequences = DATASET_OUTDIR / "sequences.db"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/complete_db.log"
    shell:
        "xonsh workflow/scripts/build_database.xsh annotate -o {output} -m {input.metadata} -h {input.chrom_names} -v {input.sv} -q {input.mc} -g {input.gffs} -c {params.column} -s {params.sp} -t {params.temp_dir} -d {params.dir} -n {params.config} {input.vcfs} -e {params.sequences} &> {log}"

        