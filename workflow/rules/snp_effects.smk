# =================================================================================================
#   Per ref_genome | Create snpeff database
# =================================================================================================


# Make symbolic links in the snpeff_data directory and create config file
rule prepare_refs_db:
    input:
        gff=REFS_DIR / "{ref_genome}" / "{ref_genome}.gff",
        fasta=rules.extract_cds_seqs.input.fasta,
        cds=rules.extract_cds_seqs.output.cds,
        prots=rules.extract_protein_seqs.output.prots,
    output:
        gff=INT_REFS_DIR / "snpeff_data" / "Species_name_{ref_genome}" / "genes.gff",
        fasta=INT_REFS_DIR / "snpeff_data" / "Species_name_{ref_genome}" / "sequences.fa",
        cds=INT_REFS_DIR / "snpeff_data" / "Species_name_{ref_genome}" / "cds.fa",
        prots=INT_REFS_DIR / "snpeff_data" / "Species_name_{ref_genome}" / "protein.fa",
    conda:
        "../envs/variants.yaml"
    params:
        name="Species_name_{ref_genome}",
        config=INT_REFS_DIR / "snpeff_data" / "snpEff.config",
    log:
        LOGS / "references" / "snpeff" / "prepare_dbs_{ref_genome}.log",
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
        gff=rules.prepare_refs_db.output.gff,
        fasta=rules.prepare_refs_db.output.fasta,
        cds=rules.prepare_refs_db.output.cds,
        prots=rules.prepare_refs_db.output.prots,
    output:
        touch(INT_REFS_DIR / "snpeff_data" / "{ref_genome}.done"),
    params:
        config=INT_REFS_DIR / "snpeff_data" / "snpEff.config",
        dir=os.getcwd() / INT_REFS_DIR / "snpeff_data",
        name="Species_name_{ref_genome}",
    log:
        LOGS / "references" / "snpeff" / "build_dbs_{ref_genome}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "snpEff build "
        "-gff3 "
        "-v "
        "-dataDir {params.dir} "
        "-config {params.config} "
        "{params.name} "
        "&>> {log}"


# =================================================================================================
#   Join samples per ref_genome | unite VCF files of all samples from each ref_genome and annotate
# =================================================================================================


rule unite_vcfs:
    input:
        unpack(unite_vcfs_input),
    output:
        vcf=INT_DATASET_DIR / "snpeff" / "{ref_genome}_union.vcf",
        tsv=INT_DATASET_DIR / "snpeff" / "{ref_genome}_presence.tsv",
    params:
        tmp_dir=os.path.join(TEMPDIR, "tmp_{ref_genome}"),
    log:
        LOGS / "dataset" / "snpeff" / "unite_vcfs_{ref_genome}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    script:
        "../scripts/unite_vcfs.xsh"


rule snpeff:
    input:
        vcf=rules.unite_vcfs.output.vcf,
        db_done=rules.build_refs_db.output,
    output:
        vcf=INT_DATASET_DIR / "snpeff" / "{ref_genome}_snpeff.vcf",
        html=INT_DATASET_DIR / "snpeff" / "{ref_genome}_snpeff.html",
    params:
        dir=os.getcwd() / INT_REFS_DIR / "snpeff_data",
        config=INT_REFS_DIR / "snpeff_data" / "snpEff.config",
        name="Species_name_{ref_genome}",
        extra=config["snpeff"]["extra"],
    log:
        LOGS / "dataset" / "snpeff" / "snpeff_{ref_genome}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "snpEff ann -v -classic "
        "-dataDir {params.dir} "
        "-config {params.config} "
        "-s {output.html} "
        "{params.extra} "
        "{params.name} "
        "{input.vcf} "
        "1> {output.vcf} 2> {log}"


rule extract_vcf_annotation:
    input:
        vcf=rules.snpeff.output.vcf,
        tsv=rules.ref_reformat_annotation.output.tsv,
        metadata=rules.quality_filter.output.metadata,
        presence=rules.unite_vcfs.output.tsv,
    output:
        effects=INT_DATASET_DIR / "snpeff" / "{ref_genome}_effects.tsv",
        variants=INT_DATASET_DIR / "snpeff" / "{ref_genome}_variants.tsv",
        lofs=INT_DATASET_DIR / "snpeff" / "{ref_genome}_lofs.tsv",
        nmds=INT_DATASET_DIR / "snpeff" / "{ref_genome}_nmds.tsv",
        classif=INT_DATASET_DIR / "snpeff" / "{ref_genome}_variant_classification.tsv",
    log:
        LOGS / "dataset" / "snpeff" / "extract_vcf_annotation_{ref_genome}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    script:
        "../scripts/extract_vcf_annotation.py"


