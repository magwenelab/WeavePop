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
LOCI_FILE = "config/loci.csv"
CHROM_NAMES = "config/chromosome_names.csv"

#### Defining variables for the reference annotation module(references.smk) ####
if config["annotate_references"]["activate"]:
    MAIN_DIR = Path(config["annotate_references"]["directory"])
    MAIN_FASTA = MAIN_DIR / config["annotate_references"]["fasta"]
    MAIN_GFF = MAIN_DIR / config["annotate_references"]["gff"]
    MAIN_NAME, _ = os.path.splitext(os.path.basename(MAIN_GFF))

#### Defining table for sample-dependent input files ####
d={'sample': SAMPLETABLE["sample"],
    'group': SAMPLETABLE["group"],
    'refgenome': REFDIR / SAMPLETABLE["group"] / (SAMPLETABLE["group"] + ".fasta"),
    'refgff': REFDIR / SAMPLETABLE["group"] / (SAMPLETABLE["group"] + ".gff")}

SAMPLE_REFERENCE = pd.DataFrame(data=d).set_index("sample", drop=False)
LINEAGE_REFERENCE = pd.DataFrame(data=d).set_index("group", drop=False)

#### Defining which final output files are being requested ####
def get_final_output():
    final_output = expand(OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa",sample=SAMPLES)
    final_output.extend(expand(OUTDIR / "snippy" / "{sample}" / "snps.bam",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}" / "lifted.gff_polished",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}" / "unmapped_features.txt",sample=SAMPLES))
    if config["annotate_references"]["activate"]:
        final_output.extend(expand(REFDIR / "{lineage}" / "{lineage}.gff",lineage=LINEAGES))
        final_output.append(REFDIR / str(MAIN_NAME + ".tsv"))
    # if config["coverage_quality"]["activate"]:
    #     final_output.extend(expand(OUTDIR / "mosdepth" / "{sample}" / "coverage.regions.bed.gz",sample=SAMPLES))
    #     final_output.extend(expand(OUTDIR / "mosdepth" / "{sample}" / "coverage_good.regions.bed.gz",sample=SAMPLES))
    #     final_output.extend(expand(OUTDIR / "mosdepth" / "{sample}" / "good_structural_variants.tsv",sample=SAMPLES))
    #     final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "distrib_cov.tsv",sample=SAMPLES))
    #     final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "mapq.bed",sample=SAMPLES))
    #     final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "mapq_window.bed",sample=SAMPLES))
    #     final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "mapq_cov_window.bed",sample=SAMPLES))
    #     final_output.extend(expand(OUTDIR / "samtools" / "{sample}" / "feature_mapq_cov.tsv",sample=SAMPLES))
    #     final_output.append(DATASET_OUTDIR / "files" / "mapping_stats.tsv")
    #     final_output.append(DATASET_OUTDIR / "files" / "structural_variants.tsv")
    #     final_output.append(DATASET_OUTDIR / "files" / "coverage_good.tsv")
    #     final_output.append(DATASET_OUTDIR / "files" / "mapqcov.tsv")
    if config["plotting"]["activate"]:
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_chrom_distribution.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_global_distribution.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_by_chrom.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_by_regions.png",sample=SAMPLES))
        final_output.append(DATASET_OUTDIR / "plots" / "dataset_depth_by_chrom.png")
        final_output.append(DATASET_OUTDIR / "plots" / "dataset_summary.png")
    if config["annotate_references"]["activate"] and config["plotting"]["activate"]:
        final_output.append(REFDIR / "unmapped_count.tsv")
        final_output.append(REFDIR / "unmapped.svg")
        final_output.append(DATASET_OUTDIR / "files" / "unmapped_count.tsv")
        final_output.append(DATASET_OUTDIR / "plots" / "unmapped.svg")
    if not config["annotate_references"]["activate"] and config["plotting"]["activate"]:
        final_output.extend(expand(DATASET_OUTDIR / "files" / "{lineage}_unmapped_count.tsv", lineage=LINEAGES))
        final_output.extend(expand(DATASET_OUTDIR / "plots" / "{lineage}_unmapped.svg", lineage=LINEAGES))
    if config["snps"]["activate"]:
        final_output.append(expand(DATASET_OUTDIR / "database.db"))
    return final_output


#### Creating softlinks to have the reference genomes in the REFDIR ####
rule links:
    input:
        REF_DATA / "{lineage}.fasta"
    output:
        REFDIR / "{lineage}" / "{lineage}.fasta"
    shell:
        "ln -s -r {input} {output}"