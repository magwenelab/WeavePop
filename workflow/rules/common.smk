# =================================================================================================
#   Load modules
# =================================================================================================

import pandas as pd
import os.path
import glob
from pathlib import Path
from snakemake.utils import validate

# =================================================================================================
#   Global variables
# =================================================================================================

UNFILTERED_SAMPLE_FILE=config["samples"]
EXCLUDE_FILE=config["samples_to_exclude"]
if EXCLUDE_FILE:
    EXCLUDE_SAMPLES=set(list(pd.read_csv(EXCLUDE_FILE, header=None, names = ["sample"])["sample"]))
    SAMPLE_ORIGINAL=(pd.read_csv(UNFILTERED_SAMPLE_FILE, sep=",", header=0))
    UNFILTERED_SAMPLE_TABLE=SAMPLE_ORIGINAL[~SAMPLE_ORIGINAL["sample"].isin(EXCLUDE_SAMPLES)]
else:
    UNFILTERED_SAMPLE_TABLE=(pd.read_csv(UNFILTERED_SAMPLE_FILE, sep=",", header=0))

validate(UNFILTERED_SAMPLE_TABLE, schema="../schemas/metadata.schema.yaml")

UNFILTERED_SAMPLES=list(set(UNFILTERED_SAMPLE_TABLE["sample"]))

REF_DATA = Path(config["references"]["directory"])
FQ_DATA = Path(config["fastqs"]["directory"])
FQ1 = config["fastqs"]["fastq_suffix1"]
FQ2 = config["fastqs"]["fastq_suffix2"]

GENERAL_OUTPUT = Path(config["output_directory"])
OUTDIR= GENERAL_OUTPUT / "samples"
DATASET_OUTDIR = GENERAL_OUTPUT / "dataset"
REFDIR = GENERAL_OUTPUT / "references"

CHROM_NAMES = config["chromosomes"]
CHROM_NAMES_TABLE = pd.read_csv(CHROM_NAMES, header=0, dtype={"chromosome": "string"})
validate(CHROM_NAMES_TABLE, schema="../schemas/chromosomes.schema.yaml")

if config["plotting"]["loci"]:
    LOCI_FILE = config["plotting"]["loci"]
    LOCI_TABLE = pd.read_csv(LOCI_FILE, sep=",", header=0)
    validate(LOCI_TABLE, schema="../schemas/loci.schema.yaml")
else:
    LOCI_FILE = OUTDIR / "loci_empty.txt"
    with open(LOCI_FILE, "w") as f:
        f.write("")
# =================================================================================================
#   Variables for the Module: Reference genomes annotation
# =================================================================================================

if config["annotate_references"]["activate"]:
    MAIN_DIR = Path(config["annotate_references"]["directory"])
    MAIN_FASTA = MAIN_DIR / config["annotate_references"]["fasta"]
    MAIN_GFF = MAIN_DIR / config["annotate_references"]["gff"]
    MAIN_NAME, _ = os.path.splitext(os.path.basename(MAIN_GFF))

# =================================================================================================
#   Tables for input functions
# =================================================================================================

d={'sample': UNFILTERED_SAMPLE_TABLE["sample"],
    'lineage': UNFILTERED_SAMPLE_TABLE["lineage"],
    'refgenome': REFDIR / UNFILTERED_SAMPLE_TABLE["lineage"] / (UNFILTERED_SAMPLE_TABLE["lineage"] + ".fasta"),
    'refgff': REFDIR / UNFILTERED_SAMPLE_TABLE["lineage"] / (UNFILTERED_SAMPLE_TABLE["lineage"] + ".gff")}

SAMPLE_REFERENCE = pd.DataFrame(data=d).set_index("sample", drop=False)
LINEAGE_REFERENCE = pd.DataFrame(data=d).set_index("lineage", drop=False)

# =================================================================================================
#   Checkpoint functions
# =================================================================================================

def listing_samples(wildcards):
    checkpoint_output = checkpoints.filtered_samples.get(**wildcards).output[0]
    return expand("{sample}", sample=glob_wildcards(os.path.join(checkpoint_output, "{sample}.txt")).sample)

SAMPLES = listing_samples

def listing_lineages(wildcards):
    checkpoint_output = checkpoints.filtered_lineages.get(**wildcards).output[0]
    return expand("{lineage}",
                lineage=glob_wildcards(os.path.join(checkpoint_output, "{lineage}.txt" )).lineage)

LINEAGES = listing_lineages

# =================================================================================================
#   Final output definition function
# =================================================================================================
def get_unfiltered_output():
    final_output = expand(OUTDIR / "snippy" / "{unf_sample}" / "snps.consensus.fa",unf_sample=UNFILTERED_SAMPLES)
    final_output.extend(expand(OUTDIR / "snippy" / "{unf_sample}" / "snps.bam",unf_sample=UNFILTERED_SAMPLES))
    final_output.extend(expand(OUTDIR / "snippy" / "{unf_sample}" / "snps.vcf.gz",unf_sample=UNFILTERED_SAMPLES))
    final_output.extend(expand(OUTDIR / "depth_quality" / "{unf_sample}" / "mapping_stats.tsv",unf_sample=UNFILTERED_SAMPLES))
    final_output.append(DATASET_OUTDIR / "depth_quality" / "mapping_stats.tsv")
    return final_output

def get_dataset_output():
    final_output = []
    final_output.append(GENERAL_OUTPUT / "metadata.csv")
    final_output.append(GENERAL_OUTPUT / "chromosomes.csv")
    final_output.append(GENERAL_OUTPUT / "loci.csv")
    final_output.append(DATASET_OUTDIR / "depth_quality" / "depth_by_chrom_good.tsv")
    final_output.append(DATASET_OUTDIR / "depth_quality" / "depth_by_chrom_raw.tsv")
    if config["cnv"]["activate"]:
        final_output.append(DATASET_OUTDIR / "cnv" / "cnv_calls.tsv")
    if config["genes_mapq_depth"]["activate"]:
        final_output.append(DATASET_OUTDIR / "depth_quality" / "feature_mapq_depth.tsv")
    if config["plotting"]["activate"]:
        final_output.append(DATASET_OUTDIR / "plots" / "dataset_depth_by_chrom.png")
        final_output.append(DATASET_OUTDIR / "plots" / "dataset_summary.png")
        if config["annotate_references"]["activate"]:
            final_output.append(REFDIR / "lifotff" / "unmapped_count.tsv")
            final_output.append(REFDIR / "lifotff" / "unmapped.svg")
            final_output.append(DATASET_OUTDIR / "liftoff" / "unmapped_count.tsv")
            final_output.append(DATASET_OUTDIR / "plots" / "unmapped.svg")
    if config["database"]["activate"]:
        final_output.append(expand(DATASET_OUTDIR / "database.db"))
    return final_output

def get_filtered_output():
    final_output = expand(OUTDIR / "liftoff" / "{sample}" / "lifted.gff_polished",sample=SAMPLES)
    final_output = final_output, expand(OUTDIR / "agat" / "{sample}" / "proteins.fa",sample=SAMPLES)
    final_output = final_output, expand(OUTDIR / "agat" / "{sample}" / "cds.fa",sample=SAMPLES)
    final_output = final_output, expand(OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv",sample=SAMPLES)
    final_output = final_output, expand(OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_raw.tsv",sample=SAMPLES)
    final_output = final_output, expand(GENERAL_OUTPUT / "filtered_lineages" / "{lineage}.txt",lineage=LINEAGES)
    if config["cnv"]["activate"]:
        final_output = final_output, expand(OUTDIR / "cnv" / "{sample}" / "cnv_calls.tsv",sample=SAMPLES)
        if config["plotting"]["activate"]:
            final_output = final_output, expand(OUTDIR / "plots" / "{sample}" / "depth_by_regions.png",sample=SAMPLES)
            final_output = final_output, expand(OUTDIR / "plots" / "{sample}" / "mapq.png",sample=SAMPLES)
    if config["genes_mapq_depth"]["activate"]:
        final_output = final_output, expand(OUTDIR / "depth_quality" / "{sample}" / "feature_mapq_depth.tsv",sample=SAMPLES)
    if config["snpeff"]["activate"]:
        final_output = final_output, expand(DATASET_OUTDIR / "snps" / "{lineage}_effects.tsv",lineage=LINEAGES)
    if config["plotting"]["activate"]:
        final_output = final_output, expand(OUTDIR / "plots" / "{sample}" / "depth_chrom_distribution.png",sample=SAMPLES)
        final_output = final_output, expand(OUTDIR / "plots" / "{sample}" / "depth_global_distribution.png",sample=SAMPLES)
        final_output = final_output, expand(OUTDIR / "plots" / "{sample}" / "depth_by_chrom.png",sample=SAMPLES)
        if not config["annotate_references"]["activate"]:
            final_output = final_output, expand(DATASET_OUTDIR / "liftoff" / "{lineage}_unmapped_count.tsv", lineage=LINEAGES)
            final_output = final_output, expand(DATASET_OUTDIR / "plots" / "{lineage}_unmapped.svg", lineage=LINEAGES)
    return final_output
    
# =================================================================================================
#   Setup rules
# =================================================================================================

# Create softlinks to have the reference genomes in the REFDIR
rule links:
    input:
        REF_DATA / "{lineage}.fasta"
    output:
        REFDIR / "{lineage}" / "{lineage}.fasta"
    shell:
        "ln -s -r {input} {output}"

rule copy_config:
    input:
        c = CHROM_NAMES,
        l = LOCI_FILE
    output:
        c = GENERAL_OUTPUT / "chromosomes.csv",
        l = GENERAL_OUTPUT / "loci.csv"
    run:
        c = pd.read_csv(input.c, header=0)
        l = pd.read_csv(input.l, header=0)
        c.to_csv(output.c, index=False)
        l.to_csv(output.l, index=False)
    
# Edit the agat config file to avoid creating log files
rule agat_config:
    output:
        REFDIR / "agat_config.yaml"
    conda:
        "../envs/agat.yaml"
    log:
        "logs/references/agat_config.log"
    shell:
        "agat config --expose &> {log} && "
        "mv agat_config.yaml {output} &> {log} && "
        "sed -i 's/log: true/log: false/g' {output} &>> {log} "

