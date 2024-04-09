
rule extract_ref_seqs:
    input:
        gff = REFDIR / "{lineage}" / "{lineage}.gff",
        fasta = REFDIR / "{lineage}" / "{lineage}.fasta"
    output:
        cds = REFDIR / "{lineage}" / "{lineage}.cds.fa",
        prots = REFDIR / "{lineage}" / "{lineage}.prots.fa"
    conda:
        "../envs/agat.yaml"
    log:
        "logs/references/extract_{lineage}.log"
    shell:
        "agat_sp_extract_sequences.pl -g {input.gff} -f {input.fasta} -o {output.cds} &> {log} "
        "&& "
        "agat_sp_extract_sequences.pl -g {input.gff} -f {input.fasta} -o {output.prots} -p &>> {log} "

rule prepare_refs_db:
    input:
        gff = rules.extract_ref_seqs.input.gff,
        fasta = rules.extract_ref_seqs.input.fasta,
        cds = rules.extract_ref_seqs.output.cds,
        prots = rules.extract_ref_seqs.output.prots
    output:
        config = DATASET_OUTDIR / "snpeff_data" / "snpEff_{lineage}.config",
        gff = DATASET_OUTDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "genes.gff",
        fasta = DATASET_OUTDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "sequences.fa",
        cds = DATASET_OUTDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "cds.fa",
        prots = DATASET_OUTDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "protein.fa"
    conda:
        "../envs/variants.yaml"
    params:
        name = config["species_name"] + "_{lineage}"
    log:
        "logs/references/prepare_dbs_{lineage}.log"
    shell:
        """
        echo "{params.name}.genome : {params.name}" > {output.config} 2> {log} && 
        ln -s -r {input.gff} {output.gff} &>> {log} && 
        ln -s -r {input.fasta} {output.fasta} &>> {log} && 
        ln -s -r {input.cds} {output.cds} &>> {log} && 
        ln -s -r {input.prots} {output.prots} &>> {log} 
        """

rule build_refs_db:
    input:
        config = rules.prepare_refs_db.output.config,
        gff = rules.prepare_refs_db.output.gff,
        fasta = rules.prepare_refs_db.output.fasta,
        cds = rules.prepare_refs_db.output.cds,
        prots = rules.prepare_refs_db.output.prots
    output:
        touch(DATASET_OUTDIR / "snpeff_data" / "{lineage}.done")
    conda:
        "../envs/variants.yaml"
    params:
        dir = os.getcwd() / DATASET_OUTDIR / "snpeff_data",
        name = config["species_name"] + "_{lineage}"
    log:
        "logs/references/build_dbs_{lineage}.log"
    shell:
        """
        snpEff build -gff3 -v -dataDir {params.dir} -config {input.config} {params.name} &>> {log}
        """


rule annotations_db:
    input:
        metadata = SAMPLEFILE,
        vcfs = expand(rules.snippy.output.vcf, sample=SAMPLES),
        databases = expand(rules.build_refs_db.output, lineage=LINEAGES)
    output:
        DATASET_OUTDIR / "annotations.db"
    params:
        column = 'group',
        sp = config["species_name"],
        temp_dir = DATASET_OUTDIR / 'tmp',
    conda:
        "../envs/variants.yaml"
    log:
        "logs/snps/annotations_db.log"
    shell:
        "xonsh workflow/scripts/build_database.xsh annotate -m {input.metadata} -c {params.column} -s {params.sp} -t {params.temp_dir} -o {output}  {input.vcfs} &> {log}"
