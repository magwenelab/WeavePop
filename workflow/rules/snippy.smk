# =================================================================================================
#   Setup rules
# =================================================================================================


rule links:
    input:
        REF_DATA / "{lineage}.fasta",
    output:
        INT_REFS_DIR / "{lineage}" / "{lineage}.fasta",
    log:
        "logs/references/links_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/shell.yaml"
    shell:
        "ln -s -r {input} {output} 2> {log}"


rule copy_chromosomes:
    input:
        CHROM_NAMES,
    output:
        DATASET_DIR / "chromosomes.csv",
    log:
        "logs/references/copy_chromosomes.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/shell.yaml"
    shell:
        """
        rsync {input} {output} 2> {log}
        """
        

# Edit the agat config file to avoid creating log files
rule agat_config:
    output:
        INT_REFS_DIR / "agat_config.yaml",
    log:
        "logs/references/agat_config.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat config --expose &> {log} && "
        "mv agat_config.yaml {output} &> {log} && "
        "sed -i 's/log: true/log: false/g' {output} &>> {log} "


# =================================================================================================
# Per sample | Run Snippy to map reads to reference genome, get assembly and call SNPs
# =================================================================================================


rule snippy:
    input:
        unpack(snippy_input),
    output:
        fa=SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.consensus.fa",
        bam=SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.bam",
        ref=SAMPLES_DIR / "snippy" / "{unf_sample}" / "ref.fa",
        bai=SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.bam.bai",
        vcf=SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.vcf.gz",
    params:
        outpath=SAMPLES_DIR / "snippy",
        tmpdir=TEMPDIR,
        extra=config["snippy"]["extra"],
    log:
        "logs/samples/snippy/snippy_{unf_sample}.log",
    threads: config["snippy"]["threads"]
    conda:
        "../envs/snippy.yaml"
    shell:
        "snippy --outdir {params.outpath}/{wildcards.unf_sample} "
        "--cpus {threads} "
        "--ref {input.refgenome} "
        "--R1 {input.fq1} "
        "--R2 {input.fq2} "
        "--force "
        "--tmpdir {params.tmpdir} "
        "{params.extra} &> {log}"
