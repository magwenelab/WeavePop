
import pandas as pd
from pathlib import Path

configfile: "config-test.yaml"


SAMPLE_TABLE = pd.read_csv(config["sample_table"]).set_index("sample", drop=False)
SAMPLES = sorted(set(SAMPLE_TABLE["sample"]))

OUTDIR = config["output_directory"]
OUTPATH = Path(OUTDIR) / "snippy"


rule snippy_all:
    input:
        expand(OUTPATH / "{sample}/snps.consensus.fa", sample=SAMPLES),
        expand(OUTPATH / "{sample}/snps.bam", sample=SAMPLES),


def snippy_from_df(wildcards):
    s = SAMPLE_TABLE.loc[wildcards.sample,]
    return {
        "fq1": s["fq1"],
        "fq2": s["fq2"],
        "refgenome": s["refgenome"],
    }


rule run_snippy:
    input:
        unpack(snippy_from_df),
    output:
        OUTPATH / "{sample}/snps.consensus.fa",
        OUTPATH / "{sample}/snps.bam",
    threads: config["threads"]
    conda:
        "envs/snippy.yaml"
    log:
        "logs/snippy/{sample}.log",
    shell:
        "snippy " 
        "--cpus {threads} "
        "--outdir {OUTPATH}/{wildcards.sample} "
        "--ref {input.refgenome} "
        "--R1 {input.fq1} "
        "--R2 {input.fq2} "
        "--force &> {log}"
