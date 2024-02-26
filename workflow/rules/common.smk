import pandas as pd
import os.path
import glob
from pathlib import Path

#### Defining global variables ####
SAMPLEFILE=config["sample_table"]
SAMPLETABLE=(pd.read_csv(config["sample_table"], sep=","))
SAMPLES=list(set(SAMPLETABLE["sample"]))
LINEAGES=list(set(SAMPLETABLE["group"]))
REF_DATA = Path(config["references"]["directory"])

FQ_DATA = Path(config["fastqs"]["directory"])
FQ1 = config["fastqs"]["fastq_suffix1"]
FQ2 = config["fastqs"]["fastq_suffix2"]

OUTDIR= Path("results/samples")
DATASET_OUTDIR = Path("results/dataset")
REFDIR = Path("results/references")

# FIX : Add conditional statements to allow for absence of these files
FEATURE_FILE = "config/features.txt"
LOCI_FILE = "config/loci.csv"
CHROM_NAMES = "config/chromosome_names.csv"

#### Defining variables for the reference annotation module(references.smk) ####
if config["annotate_references"]["activate"]:
    MAIN_DIR = Path(config["annotate_references"]["directory"])
    MAIN_FASTA = MAIN_DIR / config["annotate_references"]["fasta"]
    MAIN_GFF = MAIN_DIR / config["annotate_references"]["gff"]
    MAIN_NAME, _ = os.path.splitext(os.path.basename(MAIN_GFF))

#### Defining sample-dependent input files ####
d={'sample': SAMPLETABLE["sample"],
    'group': SAMPLETABLE["group"],
    'refgenome': REFDIR / SAMPLETABLE["group"] / (SAMPLETABLE["group"] + ".fasta"),
    'refgff': REFDIR / SAMPLETABLE["group"] / (SAMPLETABLE["group"] + ".gff")}

SAMPLE_REFERENCE = pd.DataFrame(data=d).set_index("sample", drop=False)

def snippy_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "fq1": FQ_DATA / (s["sample"] + FQ1),
        "fq2": FQ_DATA / (s["sample"] + FQ2),
        "refgenome": s["refgenome"],
    }

def liftoff_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "target": OUTDIR / "snippy" / s["sample"] / "snps.consensus.fa" ,
        "refgff": s["refgff"],
        "refgenome": s["refgenome"],
    }

def intersect_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "sampletsv": OUTDIR / "mosdepth" / s["sample"] / "ploidy_table.tsv" ,
        "maskbed": REFDIR / s["group"]  / "repeats" / "05_full" / (s["group"] + ".bed")
    }

#### Defining which final output files are being requested ####
def get_final_output():
    final_output = expand(OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa",sample=SAMPLES)
    final_output.extend(expand(OUTDIR / "snippy" / "{sample}" / "snps.bam",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}" / "lifted.gff_polished",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}" / "unmapped_features.txt",sample=SAMPLES))
    final_output.append(expand(OUTDIR / "agat" / "{sample}" / "cds.done", sample=SAMPLES))
    final_output.append(expand(OUTDIR / "agat" / "{sample}" / "prots.done", sample=SAMPLES))
    if config["annotate_references"]["activate"]:
        final_output.extend(expand(REFDIR / "{lineage}" / "{lineage}.gff",lineage=LINEAGES))
        final_output.append(REFDIR / str(MAIN_NAME + ".tsv"))
    if config["coverage_quality"]["activate"]:
        final_output.extend(expand(OUTDIR / "mosdepth" / "{sample}" / "coverage.regions.bed.gz",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "mosdepth" / "{sample}" / "coverage_good.regions.bed.gz",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "distrib_mapq.csv",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "distrib_cov.csv",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "mapq.bed",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "mapq_window.bed",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "mapq_cov_window.bed",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "annotation.gff",sample=SAMPLES))
        final_output.append(DATASET_OUTDIR / "mapping_stats.txt")
        final_output.extend(expand(OUTDIR / "mosdepth" / "{sample}" / "ploidy_table.tsv",sample=SAMPLES))
    if config["plotting"]["activate"]:
        # final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "coverage.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "cov_distribution.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "mapq_distribution.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "mapq.png",sample=SAMPLES))
        final_output.append(DATASET_OUTDIR / "plots" / "mapped_reads.png")
        final_output.append(DATASET_OUTDIR / "plots" / "cov_median_good.png")
    if config["annotate_references"]["activate"] and config["plotting"]["activate"]:
        final_output.append(REFDIR / "unmapped_count.tsv")
        final_output.append(REFDIR / "unmapped.svg")
    if not config["annotate_references"]["activate"] and config["plotting"]["activate"]:
        final_output.extend(expand(DATASET_OUTDIR / "files" / "{lineage}_unmapped_count.tsv", lineage=LINEAGES))
        final_output.extend(expand(DATASET_OUTDIR / "plots" / "{lineage}_unmapped.svg", lineage=LINEAGES))
    if config["annotate_references"]["activate"] and config["plotting"]["activate"]:
        final_output.append(DATASET_OUTDIR / "files" / "unmapped_count.tsv")
        final_output.append(DATASET_OUTDIR / "plots" / "unmapped.svg")
    
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
