
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

rule build_refs_db:
    input:
        rules.extract_ref_seqs.input.gff,
        rules.extract_ref_seqs.input.fasta,
        rules.extract_ref_seqs.output.cds,
        rules.extract_ref_seqs.output.prots
    output:
        touch(DATASET_OUTDIR / "snpeff_{lineage}.done")
    conda:
        "../envs/variants.yaml"
    params:
        config["species_name"] + "_{lineage}"
    log:
        "logs/references/build_dbs_{lineage}.log"
    shell:
        """
        echo {params}.genome : {params} >> $CONDA_PREFIX/share/snpeff-5.2-0/snpEff.config 2> {log}
        mkdir -p $CONDA_PREFIX/share/snpeff-5.2-0/data/{params} &>> {log}
        scp {input[0]} $CONDA_PREFIX/share/snpeff-5.2-0/data/{params}/genes.gff &>> {log}
        scp {input[1]} $CONDA_PREFIX/share/snpeff-5.2-0/data/{params}/sequences.fa &>> {log}
        scp {input[2]} $CONDA_PREFIX/share/snpeff-5.2-0/data/{params}/cds.fa &>> {log}
        scp {input[3]} $CONDA_PREFIX/share/snpeff-5.2-0/data/{params}/protein.fa &>> {log}
        snpEff build -gff3 -v {params} &>> {log}
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
