# =================================================================================================
#   Load modules
# =================================================================================================

import pandas as pd
import os.path
import glob
from pathlib import Path
from snakemake.utils import validate

# =================================================================================================
#   Define global variables and validate input files using schemas
# =================================================================================================

UNFILT_SAMPLE_FILE = config["metadata"]
EXCLUDE_FILE = config["samples_to_exclude"]
if EXCLUDE_FILE:
    EXCLUDE_SAMPLES = set(list(pd.read_csv(EXCLUDE_FILE, header=None, names=["sample"])["sample"]))
    SAMPLE_ORIGINAL = pd.read_csv(UNFILT_SAMPLE_FILE, sep=",", header=0)
    UNFILT_SAMPLE_TABLE = SAMPLE_ORIGINAL[~SAMPLE_ORIGINAL["sample"].isin(EXCLUDE_SAMPLES)]
else:
    UNFILT_SAMPLE_TABLE = pd.read_csv(UNFILT_SAMPLE_FILE, sep=",", header=0)

validate(UNFILT_SAMPLE_TABLE, schema="../schemas/metadata.schema.yaml")

UNFILT_SAMPLES = list(set(UNFILT_SAMPLE_TABLE["sample"]))

REF_DATA = Path(config["references"]["directory"])
FQ_DATA = Path(config["fastqs"]["directory"])
FQ1 = config["fastqs"]["fastq_suffix1"]
FQ2 = config["fastqs"]["fastq_suffix2"]

OUTPUT = Path(config["output_directory"])
SAMPLES_DIR = OUTPUT / SAMPLES_DIR_NAME
DATASET_DIR = OUTPUT / DATASET_DIR_NAME
REFS_DIR = OUTPUT / REFS_DIR_NAME
INTDIR = OUTPUT / INTDIR_NAME
INT_SAMPLES_DIR = INTDIR / SAMPLES_DIR_NAME
INT_DATASET_DIR = INTDIR / DATASET_DIR_NAME
INT_REFS_DIR = INTDIR / REFS_DIR_NAME
TEMPDIR = str(INTDIR / TEMPDIR_NAME)

CHROM_NAMES = config["chromosomes"]
CHROM_NAMES_TABLE = pd.read_csv(CHROM_NAMES, header=0, dtype={"chromosome": "string"})
validate(CHROM_NAMES_TABLE, schema="../schemas/chromosomes.schema.yaml")

if config["plotting"]["loci"]:
    LOCI_FILE = config["plotting"]["loci"]
    LOCI_TABLE = pd.read_csv(LOCI_FILE, sep=",", header=0)
    validate(LOCI_TABLE, schema="../schemas/loci.schema.yaml")
else:
    LOCI_FILE = SAMPLES_DIR / "loci_empty.txt"
    with open(LOCI_FILE, "w") as f:
        f.write("")

# =================================================================================================
#   Variables for the module Reference annotation
# =================================================================================================

if config["annotate_references"]["activate"]:
    MAIN_DIR = Path(config["annotate_references"]["directory"])
    MAIN_FASTA = MAIN_DIR / config["annotate_references"]["fasta"]
    MAIN_GFF = MAIN_DIR / config["annotate_references"]["gff"]
    MAIN_NAME, _ = os.path.splitext(os.path.basename(MAIN_GFF))

# =================================================================================================
#   Tables for input functions
# =================================================================================================

d = {
    "sample": UNFILT_SAMPLE_TABLE["sample"],
    "lineage": UNFILT_SAMPLE_TABLE["lineage"],
    "refgenome": INT_REFS_DIR
    / UNFILT_SAMPLE_TABLE["lineage"]
    / (UNFILT_SAMPLE_TABLE["lineage"] + ".fasta"),
    "refgff": INT_REFS_DIR / UNFILT_SAMPLE_TABLE["lineage"] 
    / (UNFILT_SAMPLE_TABLE["lineage"] + "_interg_introns.gff"),
}

SAMPLE_REFERENCE = pd.DataFrame(data=d).set_index("sample", drop=False)
LINEAGE_REFERENCE = pd.DataFrame(data=d).set_index("lineage", drop=False)

# =================================================================================================
#   Input functions for rules
# =================================================================================================

def snippy_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.unf_sample,]
    return {
        "fq1": FQ_DATA / (s["sample"] + FQ1),
        "fq2": FQ_DATA / (s["sample"] + FQ2),
        "refgenome": s["refgenome"],
    }


def liftoff_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "target": SAMPLES_DIR / "snippy" / s["sample"] / "snps.consensus.fa",
        "refgff": s["refgff"],
        "refgenome": s["refgenome"],
    }


def depth_distribution_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.unf_sample,]
    return {
        "bam": SAMPLES_DIR / "snippy" / s["sample"] / "snps.bam",
        "bai": SAMPLES_DIR / "snippy" / s["sample"] / "snps.bam.bai",
        "bam_good": INT_SAMPLES_DIR / "depth_quality" / s["sample"] / "snps_good.bam",
        "bai_good": INT_SAMPLES_DIR / "depth_quality" / s["sample"] / "snps_good.bam.bai",
    }


def depth_by_windows_plots_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "depth": INT_SAMPLES_DIR / "depth_quality" / s["sample"] / "depth_by_windows.tsv",
        "cnv": SAMPLES_DIR / "cnv" / s["sample"] / "cnv_calls.tsv",
        "repeats": REFS_DIR / s["lineage"] / (s["lineage"] + "_repeats.bed"),
    }


def mapq_plot_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "mapq": INT_SAMPLES_DIR / "depth_quality" / s["sample"] / "mapq_window.bed",
        "cnv": SAMPLES_DIR / "cnv" / s["sample"] / "cnv_calls.tsv",
        "repeats": REFS_DIR / s["lineage"] / (s["lineage"] + "_repeats.bed"),
    }


def cnv_calling_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "depth": INT_SAMPLES_DIR / "depth_quality" / s["sample"] / "depth_by_windows.tsv",
        "repeats": REFS_DIR / s["lineage"] / (s["lineage"] + "_repeats.bed"),
    }


def intersect_vcfs_input(wildcards):
    sample_wildcards = listing_samples(wildcards)
    l = LINEAGE_REFERENCE[LINEAGE_REFERENCE["sample"].isin(sample_wildcards)]
    l = l.loc[wildcards.lineage,]
    return {"vcfs": expand(SAMPLES_DIR / "snippy" / "{sample}" / "snps.vcf.gz", sample=l["sample"])}


# =================================================================================================
#   Checkpoint functions
# =================================================================================================

def listing_samples(wildcards):
    checkpoint_output = checkpoints.filter_wildcards.get(**wildcards).output[0]
    return expand(
        "{sample}", sample=glob_wildcards(os.path.join(checkpoint_output, "{sample}.txt")).sample
    )

SAMPLES = listing_samples


def listing_lineages(wildcards):
    checkpoint_output = checkpoints.filter_wildcards.get(**wildcards).output[1]
    return expand(
        "{lineage}",
        lineage=glob_wildcards(os.path.join(checkpoint_output, "{lineage}.txt")).lineage,
    )

LINEAGES = listing_lineages


# =================================================================================================
#   Final output definition functions
# =================================================================================================

# --Output per sample previous to sample filtering-------------------------------------------------
def get_unfiltered_output():
    final_output = expand(
        SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.consensus.fa", unf_sample=UNFILT_SAMPLES
    )
    final_output.extend(
        expand(SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.bam", unf_sample=UNFILT_SAMPLES)
    )
    final_output.extend(
        expand(SAMPLES_DIR / "snippy" / "{unf_sample}" / "snps.vcf.gz", unf_sample=UNFILT_SAMPLES)
    )
    return final_output


# --Output per sample after sample filtering-------------------------------------------------------
def get_filtered_output():
    final_output = expand(
        SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff", sample=SAMPLES
    )
    final_output = final_output, expand(
        SAMPLES_DIR / "annotation" / "{sample}" / "proteins.fa", sample=SAMPLES
    )
    final_output = final_output, expand(
        SAMPLES_DIR / "annotation" / "{sample}" / "cds.fa", sample=SAMPLES
    )
    final_output = final_output, expand(
        INT_REFS_DIR / "filtered_lineages" / "{lineage}.txt", lineage=LINEAGES
    )
    if config["plotting"]["activate"]:
        final_output = final_output, expand(
            SAMPLES_DIR / "plots" / "{sample}" / "depth_by_windows.png", sample=SAMPLES
        )
        final_output = final_output, expand(
            SAMPLES_DIR / "plots" / "{sample}" / "mapq.png", sample=SAMPLES
        )
        final_output = final_output, expand(
            SAMPLES_DIR / "plots" / "{sample}" / "depth_chrom_distribution.png", sample=SAMPLES
        )
        final_output = final_output, expand(
            SAMPLES_DIR / "plots" / "{sample}" / "depth_global_distribution.png", sample=SAMPLES
        )
        final_output = final_output, expand(
            SAMPLES_DIR / "plots" / "{sample}" / "depth_by_chrom.png", sample=SAMPLES
        )
    return final_output


# --Output per dataset---------------------------------------------------------------------------
def get_dataset_output():
    final_output = []
    final_output.append(DATASET_DIR / "metadata.csv")
    final_output.append(DATASET_DIR / "chromosomes.csv")
    final_output.append(DATASET_DIR / "depth_quality" / "mapping_stats.tsv")
    if config["annotate_references"]["activate"]:
        final_output.append(REFS_DIR / "refs_unmapped_features.tsv")
    if config["genes_mapq_depth"]["activate"]:
        final_output.append(DATASET_DIR / "depth_quality" / "feature_mapq_depth.tsv")
    if config["snpeff"]["activate"]:
        final_output.append(DATASET_DIR / "snps" / "effects.tsv")
    if config["cnv"]["activate"]:
        final_output.append(DATASET_DIR / "cnv" / "cnv_calls.tsv")
    if config["plotting"]["activate"]:
        final_output.append(DATASET_DIR / "plots" / "dataset_depth_by_chrom.png")
        final_output.append(DATASET_DIR / "plots" / "dataset_summary.png")
    if config["database"]["activate"]:
        final_output.append(expand(DATASET_DIR / "database.db"))
    return final_output
