
import pandas as pd
from pathlib import Path


configfile: "config-agat.yaml"


SAMPLE_TABLE = pd.read_csv(config["sample_table"]).set_index("sample", drop=False)
SAMPLES = list(set(SAMPLE_TABLE["sample"]))

OUTDIR = config["output_directory"]
OUTPATH = Path(OUTDIR) / "agat"



rule agat_all:
    input:
        expand(OUTPATH / "{sample}/predicted_cds.fa", sample=SAMPLES),
        expand(OUTPATH / "{sample}/predicted_proteins.fa", sample=SAMPLES),


def agat_from_df(wildcards):
    s = SAMPLE_TABLE.loc[wildcards.sample,]
    return {
        "gff": s["samplegff"],
        "fasta": s["samplegenome"],
    }


rule run_agat:
    input:
        unpack(agat_from_df)
    output:
        cds = OUTPATH / "{sample}/predicted_cds.fa",
        prots = OUTPATH / "{sample}/predicted_proteins.fa"
    conda:
        "envs/agat.yaml"
    log: 
        cds = "logs/agat/{sample}_cds.log",
        prots = "logs/agat/{sample}_prots.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fasta} "
        "-o {output.cds} "
        "&> {log.cds} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.prots} "
        "-p  &> {log.prots} " 
        " && "
        "rm lifted.agat.log || true"