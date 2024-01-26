import pandas as pd
import os.path
import glob
from pathlib import Path

#### Defining global variables ####
sampletable=(pd.read_csv(config["sample_table"], sep=","))
SAMPLES=list(set(sampletable["sample"]))
LINEAGES=list(set(sampletable["group"]))
REF_DATA = Path(config["references"]["directory"])

FQ_DATA = Path(config["fastqs"]["directory"])
FQ1 = config["fastqs"]["fastq_suffix1"]
FQ2 = config["fastqs"]["fastq_suffix2"]

OUTDIR= Path("results/samples")
DATASET_OUTDIR = Path("results/dataset")
REFDIR = Path("results/references")

FEATURE_FILE = "config/features.txt"
#### Defining variables for the reference annotation module(references.smk) ####
if config["annotate_references"]["activate"]:
    MAIN_DIR = Path(config["annotate_references"]["directory"])
    MAIN_FASTA = MAIN_DIR / config["annotate_references"]["fasta"]
    MAIN_GFF = MAIN_DIR / config["annotate_references"]["gff"]
    MAIN_NAME, _ = os.path.splitext(os.path.basename(MAIN_GFF))

#### Defining sample-dependent input files ####
d={'sample': sampletable["sample"],
    'group': sampletable["group"],
    'fq1': FQ_DATA / (sampletable["sample"] + FQ1),
    'fq2': FQ_DATA / (sampletable["sample"] + FQ2),
    'refgenome' : REFDIR / sampletable["group"] / (sampletable["group"] + ".fasta"),
    'refgff' : REFDIR / sampletable["group"] / (sampletable["group"] + ".gff")}
SAMPLE_REFERENCE = pd.DataFrame(data=d).set_index("sample", drop=False)
def snippy_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "fq1": s["fq1"],
        "fq2": s["fq2"],
        "refgenome": s["refgenome"],
    }

def liftoff_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "target": OUTDIR / "snippy" / s["sample"] / "snps.consensus.fa" ,
        "refgff": s["refgff"],
        "refgenome": s["refgenome"],
    }

#### Defining which final output files are being requested ####
def get_final_output():
    final_output = expand(OUTDIR / "snippy" / "{sample}/snps.consensus.fa",sample=SAMPLES)
    final_output.extend(expand(OUTDIR / "snippy" / "{sample}/snps.bam",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}/lifted.gff_polished",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}/unmapped_features.txt",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "agat" / "{sample}/cds.fa",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "agat" / "{sample}/proteins.fa",sample=SAMPLES))
    final_output.append(DATASET_OUTDIR / "cds.done")
    final_output.append(DATASET_OUTDIR / "prots.done")
    if config["annotate_references"]["activate"]:
        final_output.extend(expand(REFDIR / "{lineage}" / "{lineage}.gff",lineage=LINEAGES))
        final_output.append(REFDIR / str(MAIN_NAME + ".tsv"))
    return final_output

#### Creating softlinks to have the reference genomes in the REFDIR ####
rule links:
    input:
        REF_DATA / "{lineage}.fasta"
    output:
        REFDIR / "{lineage}" / "{lineage}.fasta"
    shell:
        "ln -s -r {input} {output}"

if not config["annotate_references"]["activate"]:
    rule gff_links:
        input:
            REF_DATA / "{lineage}.gff"
        output:
            REFDIR / "{lineage}" / "{lineage}.gff"
        shell:
            "ln -s -r {input} {output}"
