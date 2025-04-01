# =================================================================================================
#   Per lineage | Create snpeff database
# =================================================================================================

# Make symbolic links in the snpeff_data directory and create config file
rule prepare_refs_db:
    input:
        gff=INT_REFS_DIR / "{lineage}" / "{lineage}_interg_introns.gff",
        fasta=rules.extract_cds_seqs.input.fasta,
        cds=rules.extract_cds_seqs.output.cds,
        prots=rules.extract_protein_seqs.output.prots,
    output:
        gff=INT_REFS_DIR / "snpeff_data" / "Species_name_{lineage}" / "genes.gff",
        fasta=INT_REFS_DIR
        / "snpeff_data"
        / "Species_name_{lineage}"
        / "sequences.fa",
        cds=INT_REFS_DIR / "snpeff_data" / "Species_name_{lineage}" / "cds.fa",
        prots=INT_REFS_DIR
        / "snpeff_data"
        / "Species_name_{lineage}"
        / "protein.fa",
    conda:
        "../envs/variants.yaml"
    params:
        name="Species_name_{lineage}",
        config=INT_REFS_DIR / "snpeff_data" / "snpEff.config",
    log:
        LOGS / "references" / "snpeff" / "prepare_dbs_{lineage}.log",
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
        name="Species_name_{lineage}",
    log:
        LOGS / "references" / "snpeff" / "build_dbs_{lineage}.log",
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
        vcf=INT_DATASET_DIR / "snpeff" / "{lineage}_intersection.vcf",
        tsv=INT_DATASET_DIR / "snpeff" / "{lineage}_presence.tsv",
    params:
        tmp_dir=os.path.join(TEMPDIR, "tmp_{lineage}"),
    log:
        LOGS / "dataset" / "snpeff" / "intersect_vcfs_{lineage}.log",
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
        vcf=INT_DATASET_DIR / "snpeff" / "{lineage}_snpeff.vcf",
        html=INT_DATASET_DIR / "snpeff" / "{lineage}_snpeff.html",
    params:
        dir=os.getcwd() / INT_REFS_DIR / "snpeff_data",
        config=INT_REFS_DIR / "snpeff_data" / "snpEff.config",
        name="Species_name_{lineage}",
    log:
        LOGS / "dataset" / "snpeff" / "snpeff_{lineage}.log",
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
        tsv=rules.ref_gff2tsv.output.tsv,
    output:
        effects=INT_DATASET_DIR / "snpeff" / "{lineage}_effects.tsv",
        variants=INT_DATASET_DIR / "snpeff" / "{lineage}_variants.tsv",
        lofs=INT_DATASET_DIR / "snpeff" / "{lineage}_lofs.tsv",
        nmds=INT_DATASET_DIR / "snpeff" / "{lineage}_nmds.tsv",
    log:
        LOGS / "dataset" / "snpeff" / "extract_vcf_annotation_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "xonsh workflow/scripts/extract_vcf_annotation.xsh "
        "-i {input.vcf} "
        "-g {input.tsv} "
        "-e {output.effects} "
        "-v {output.variants} "
        "-f {output.lofs} "
        "-n {output.nmds} "
        "-l {wildcards.lineage} "
        "&> {log}"
