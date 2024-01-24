import pandas as pd
from pathlib import Path


configfile: "config-test.yaml"


SAMPLE_TABLE = pd.read_csv(config["sample_table"]).set_index("sample", drop=False)
SAMPLES = list(set(SAMPLE_TABLE["sample"]))


REFDIR = config["reference_directory"]
OUTDIR = config["output_directory"]
FASTQDIR = config["fastq_directory"]

REFPATH, OUTPATH, FASTQPATH = [
    Path(config[i])
    for i in ("reference_directory", "output_directory", "fastq_directory")
]


rule:
    input:
        expand(OUTPATH / "{sample}/snps.consensus.fa", sample=SAMPLES),
        expand(OUTPATH / "{sample}/snps.bam", sample=SAMPLES),


def snippy_from_df(wildcards):
    s = SAMPLE_TABLE.loc[wildcards.sample,]
    return {
        "fq1": FASTQPATH / s["fq1"],
        "fq2": FASTQPATH / s["fq2"],
        "refgenome": REFPATH / s["refgenome"],
    }


rule snippy:
    input:
        unpack(snippy_from_df),
    output:
        OUTPATH / "{sample}/snps.consensus.fa",
        OUTPATH / "{sample}/snps.bam",
    threads: config["threads_snippy"]
    conda:
        "envs/snippy.yaml"
    log:
        "logs/snippy/{sample}.log",
    shell:
        "snippy --outdir {OUTPATH}/{wildcards.sample} "
        "--cpus {threads} "
        "--ref {input.refgenome} "
        "--R1 {input.fq1} "
        "--R2 {input.fq2} "
        "--force &> {log}"
