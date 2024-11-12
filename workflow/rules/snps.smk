# =================================================================================================
#   Per lineage | Extract CDS and protein sequences from reference genomes
# =================================================================================================
# Extract cds and protein sequences from reference genomes
rule extract_cds_seqs:
    input:
        gff=REFS_DIR / "{lineage}" / "{lineage}.gff",
        fasta=INT_REFS_DIR / "{lineage}" / "{lineage}.fasta",
        config=rules.agat_config.output,
    output:
        cds=INT_REFS_DIR / "{lineage}" / "{lineage}.cds.fa",
    log:
        "logs/references/extract_cds_seqs_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} "
        "-f {input.fasta} "
        "-o {output.cds} "
        "-c {input.config} "
        "&> {log}"


rule extract_protein_seqs:
    input:
        gff=REFS_DIR / "{lineage}" / "{lineage}.gff",
        fasta=INT_REFS_DIR / "{lineage}" / "{lineage}.fasta",
        config=rules.agat_config.output,
        cds=rules.extract_cds_seqs.output.cds,
    output:
        prots=INT_REFS_DIR / "{lineage}" / "{lineage}.prots.fa",
    log:
        "logs/references/extract_protein_seqs_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} "
        "-f {input.fasta} "
        "-o {output.prots} "
        "-p "
        "-c {input.config} "
        "&> {log}"


# Make symbolic links in the snpeff_data directory and create config file
rule prepare_refs_db:
    input:
        gff=rules.extract_cds_seqs.input.gff,
        fasta=rules.extract_cds_seqs.input.fasta,
        cds=rules.extract_cds_seqs.output.cds,
        prots=rules.extract_protein_seqs.output.prots,
    output:
        gff=INT_REFS_DIR
        / "snpeff_data"
        / str(config["species_name"] + "_{lineage}")
        / "genes.gff",
        fasta=INT_REFS_DIR
        / "snpeff_data"
        / str(config["species_name"] + "_{lineage}")
        / "sequences.fa",
        cds=INT_REFS_DIR
        / "snpeff_data"
        / str(config["species_name"] + "_{lineage}")
        / "cds.fa",
        prots=INT_REFS_DIR
        / "snpeff_data"
        / str(config["species_name"] + "_{lineage}")
        / "protein.fa",
    conda:
        "../envs/variants.yaml"
    params:
        name=config["species_name"] + "_{lineage}",
        config=INT_REFS_DIR / "snpeff_data" / "snpEff.config",
    log:
        "logs/references/prepare_dbs_{lineage}.log",
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
        touch(INT_REFS_DIR / "snpeff_data" / "{lineage}.done"),
    params:
        config=INT_REFS_DIR / "snpeff_data" / "snpEff.config",
        dir=os.getcwd() / INT_REFS_DIR / "snpeff_data",
        name=config["species_name"] + "_{lineage}",
    log:
        "logs/references/build_dbs_{lineage}.log",
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
#   Join samples per lineage | Intersect VCF files of all samples from each lineage and annotate
# =================================================================================================


rule intersect_vcfs:
    input:
        unpack(intersect_vcfs_input),
    output:
        vcf=INT_DATASET_DIR / "snps" / "{lineage}_intersection.vcf",
        tsv=INT_DATASET_DIR / "snps" / "{lineage}_presence.tsv",
    params:
        tmp_dir=os.path.join(TEMPDIR, "tmp_{lineage}"),
    log:
        "logs/dataset/snps/intersect_vcfs_{lineage}.log",
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
        db_done=rules.build_refs_db.output,
    output:
        vcf=INT_DATASET_DIR / "snps" / "{lineage}_snpeff.vcf",
        html=INT_DATASET_DIR / "snps" / "{lineage}_snpeff.html",
    params:
        dir=os.getcwd() / INT_REFS_DIR / "snpeff_data",
        config=INT_REFS_DIR / "snpeff_data" / "snpEff.config",
        name=config["species_name"] + "_{lineage}",
    log:
        "logs/dataset/snps/snpeff_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "snpEff ann -v -classic "
        "-dataDir {params.dir} "
        "-config {params.config} "
        "-s {output.html} "
        "{params.name} "
        "{input.vcf} "
        "1> {output.vcf} 2> {log}"


rule extract_vcf_annotation:
    input:
        vcf=rules.snpeff.output.vcf,
        gff=rules.ref_reformat_annotation.output.tsv,
    output:
        effects=INT_DATASET_DIR / "snps" / "{lineage}_effects.tsv",
        variants=INT_DATASET_DIR / "snps" / "{lineage}_variants.tsv",
        lofs=INT_DATASET_DIR / "snps" / "{lineage}_lofs.tsv",
        nmds=INT_DATASET_DIR / "snps" / "{lineage}_nmds.tsv",
    log:
        "logs/dataset/snps/extract_vcf_annotation_{lineage}.log",
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
