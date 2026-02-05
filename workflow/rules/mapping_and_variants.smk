# =================================================================================================
# Per sample | Map reads to reference genome, get assembly and call SNPs
# =================================================================================================


rule mapping_and_variants:
    input:
        unpack(mapping_and_variants_input),
    output:
        fa=SAMPLES_DIR / "mapping_and_variants" / "{unf_sample}" / "snps.consensus.fa",
        bam=SAMPLES_DIR / "mapping_and_variants" / "{unf_sample}" / "snps.bam",
        ref=SAMPLES_DIR / "mapping_and_variants" / "{unf_sample}" / "ref.fa",
        bai=SAMPLES_DIR / "mapping_and_variants" / "{unf_sample}" / "snps.bam.bai",
        vcf=SAMPLES_DIR / "mapping_and_variants" / "{unf_sample}" / "snps.vcf.gz",
    params:
        outpath=SAMPLES_DIR / "mapping_and_variants",
        tmpdir=TEMPDIR,
        extra=config["mapping_and_variants"]["extra"],
    log:
        LOGS
        / "samples"
        / "mapping_and_variants"
        / "mapping_and_variants_{unf_sample}.log",
    threads: config["mapping_and_variants"]["threads"]
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
