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

SAMPLEFILE=config["samples"]
SAMPLETABLE=(pd.read_csv(config["samples"], sep=",", header=0))
validate(SAMPLETABLE, schema="../schemas/metadata.schema.yaml")

SAMPLES=list(set(SAMPLETABLE["sample"]))
LINEAGES=list(set(SAMPLETABLE["lineage"]))
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
#   Tables for sample-dependent input files 
# =================================================================================================

d={'sample': SAMPLETABLE["sample"],
    'lineage': SAMPLETABLE["lineage"],
    'refgenome': REFDIR / SAMPLETABLE["lineage"] / (SAMPLETABLE["lineage"] + ".fasta"),
    'refgff': REFDIR / SAMPLETABLE["lineage"] / (SAMPLETABLE["lineage"] + ".gff")}

SAMPLE_REFERENCE = pd.DataFrame(data=d).set_index("sample", drop=False)
LINEAGE_REFERENCE = pd.DataFrame(data=d).set_index("lineage", drop=False)

# =================================================================================================
#   Final output definition function
# =================================================================================================

def get_final_output():
    final_output = expand(OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa",sample=SAMPLES)
    final_output.extend(expand(OUTDIR / "snippy" / "{sample}" / "snps.bam",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}" / "lifted.gff_polished",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}" / "unmapped_features.txt",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "agat" / "{sample}" / "proteins.fa",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "agat" / "{sample}" / "cds.fa",sample=SAMPLES))
    final_output.append(GENERAL_OUTPUT / "metadata.csv")
    final_output.append(GENERAL_OUTPUT / "chromosomes.csv")
    final_output.append(GENERAL_OUTPUT / "loci.csv")
    if config["depth_quality"]["activate"]:
        final_output.extend(expand(OUTDIR / "depth_quality" / "{sample}" / "depth_by_regions.tsv",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_raw.tsv",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "depth_quality" / "{sample}" / "mapping_stats.tsv",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "depth_quality" / "{sample}" / "feature_mapq_depth.tsv",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "cnv" / "{sample}" / "cnv_calls.tsv",sample=SAMPLES))
    if config["snpeff"]["activate"]:
        final_output.extend(expand(DATASET_OUTDIR / "snps" / "{lineage}_effects.tsv",lineage=LINEAGES))
    if config["plotting"]["activate"]:
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_chrom_distribution.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_global_distribution.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_by_chrom.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_by_regions.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "mapq.png",sample=SAMPLES))
    if config["annotate_references"]["activate"] and config["plotting"]["activate"]:
        final_output.append(REFDIR / "lifotff" / "unmapped_count.tsv")
        final_output.append(REFDIR / "lifotff" / "unmapped.svg")
        final_output.append(DATASET_OUTDIR / "liftoff" / "unmapped_count.tsv")
        final_output.append(DATASET_OUTDIR / "plots" / "unmapped.svg")
    if not config["annotate_references"]["activate"] and config["plotting"]["activate"]:
        final_output.extend(expand(DATASET_OUTDIR / "liftoff" / "{lineage}_unmapped_count.tsv", lineage=LINEAGES))
        final_output.extend(expand(DATASET_OUTDIR / "plots" / "{lineage}_unmapped.svg", lineage=LINEAGES))
    if config["database"]["activate"]:
        final_output.append(DATASET_OUTDIR / "depth_quality" / "feature_mapq_depth.tsv")
        final_output.append(DATASET_OUTDIR / "cnv" / "cnv_calls.tsv")
        final_output.append(expand(DATASET_OUTDIR / "database.db"))
        if config["plotting"]["activate"]:
            final_output.append(DATASET_OUTDIR / "plots" / "dataset_depth_by_chrom.png")
            final_output.append(DATASET_OUTDIR / "plots" / "dataset_summary.png")
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
        s = config["samples"],
        c = config["chromosomes"],
        l = LOCI_FILE
    output:
        s = GENERAL_OUTPUT / "metadata.csv",
        c = GENERAL_OUTPUT / "chromosomes.csv",
        l = GENERAL_OUTPUT / "loci.csv"
    shell:
        "scp {input.s} {output.s} && "
        "scp {input.c} {output.c} && "
        "scp {input.l} {output.l}"
    
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

