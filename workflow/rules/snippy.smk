# =================================================================================================
# Per sample | Run Snippy to map reads to reference genome, get assembly and call SNPs
# =================================================================================================

def snippy_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.unf_sample,]
    return {
        "fq1": FQ_DATA / (s["sample"] + FQ1),
        "fq2": FQ_DATA / (s["sample"] + FQ2),
        "refgenome": s["refgenome"],
    }

rule snippy:
    input:
        unpack(snippy_input)
    output:
        fa = SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.consensus.fa",
        bam = SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.bam",
        ref = SAMPLES_DIR / "snippy" / "{unf_sample}" / "ref.fa",
        bai = SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.bam.bai",
        vcf = SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.vcf.gz"
    threads: 
        config["snippy"]["threads"]
    params:
        outpath = SAMPLES_DIR / "snippy",
        extra = config["snippy"]["extra"]
    conda:
        "../envs/snippy.yaml"
    log:
        "logs/samples/snippy/snippy_{unf_sample}.log"
    shell:
        "snippy --outdir {params.outpath}/{wildcards.unf_sample} "
        "--cpus {threads} "
        "--ref {input.refgenome} "
        "--R1 {input.fq1} "
        "--R2 {input.fq2} "
        "--force "
        "{params.extra} &> {log}"
