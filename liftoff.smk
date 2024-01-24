
import pandas as pd
from pathlib import Path


configfile: "config-liftoff.yaml"


SAMPLE_TABLE = pd.read_csv(config["sample_table"]).set_index("sample", drop=False)
SAMPLES = list(set(SAMPLE_TABLE["sample"]))

OUTDIR = config["output_directory"]
OUTPATH = Path(OUTDIR) / "liftoff"


FEATURE_FILE = config["feature_file"]


rule liftoff_all:
    input:
        expand(OUTPATH / "{sample}/lifted.gff_polished", sample=SAMPLES),
        expand(OUTPATH / "{sample}/snps.bam", sample=SAMPLES),


def liftoff_from_df(wildcards):
    s = SAMPLE_TABLE.loc[wildcards.sample,]
    return {
        "samplegenome": s["samplegenome"],
        "refgff": s["refgff"],
        "refgenome": s["refgenome"],
    }


rule run_liftoff:
    input:
        unpack(liftoff_from_df),
        features=FEATURE_FILE,
    output:
        OUTPATH / "{sample}/lifted.gff",
        OUTPATH / "{sample}/lifted.gff_polished",
        OUTPATH / "{sample}/unmapped_features.txt",
    threads: config["threads"]
    conda:
        "envs/liftoff.yaml"
    log:
        "logs/liftoff/{sample}.log",
    shell:
        "liftoff "
        "-g {input.refgff} "
        "-polish "
        "-f {input.features} "
        "-dir {OUTPATH}/{wildcards.sample}/intermediate_files "
        "-u {OUTPATH}/{wildcards.sample}/unmapped_features.txt "
        "-o {OUTPATH}/{wildcards.sample}/lifted.gff "
        "-p {threads} "
        "{input.samplegenome} "
        "{input.refgenome} &> {log}"
