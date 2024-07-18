# =================================================================================================
#   Load modules
# =================================================================================================

import pandas as pd
import os.path
import glob
from pathlib import Path

# =================================================================================================
#   Global variables
# =================================================================================================

SAMPLEFILE=config["sample_table"]
SAMPLETABLE=(pd.read_csv(config["sample_table"], sep=","))
SAMPLES=list(set(SAMPLETABLE["sample"]))
LINEAGES=list(set(SAMPLETABLE["group"]))
REF_DATA = Path(config["references"]["directory"])

FQ_DATA = Path(config["fastqs"]["directory"])
FQ1 = config["fastqs"]["fastq_suffix1"]
FQ2 = config["fastqs"]["fastq_suffix2"]

GENERAL_OUTPUT = Path(config["output_directory"])
OUTDIR= GENERAL_OUTPUT / "samples"
DATASET_OUTDIR = GENERAL_OUTPUT / "dataset"
REFDIR = GENERAL_OUTPUT / "references"

CHROM_NAMES = config["chromosome_names"]
LOCI_FILE = config["loci"]

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
    'group': SAMPLETABLE["group"],
    'refgenome': REFDIR / SAMPLETABLE["group"] / (SAMPLETABLE["group"] + ".fasta"),
    'refgff': REFDIR / SAMPLETABLE["group"] / (SAMPLETABLE["group"] + ".gff")}

SAMPLE_REFERENCE = pd.DataFrame(data=d).set_index("sample", drop=False)
LINEAGE_REFERENCE = pd.DataFrame(data=d).set_index("group", drop=False)

# =================================================================================================
#   Final output definition function
# =================================================================================================

def get_final_output():
    final_output = expand(OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa",sample=SAMPLES)
    final_output.extend(expand(OUTDIR / "snippy" / "{sample}" / "snps.bam",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}" / "lifted.gff_polished",sample=SAMPLES))
    final_output.extend(expand(OUTDIR / "liftoff" / "{sample}" / "unmapped_features.txt",sample=SAMPLES))
    if config["coverage_quality"]["activate"]:
        final_output.append(DATASET_OUTDIR / "files" / "mapping_stats.tsv")
        final_output.append(DATASET_OUTDIR / "files" / "cnv_calls.tsv")
        final_output.append(DATASET_OUTDIR / "files" / "depth_by_chrom_good.tsv")
        final_output.append(DATASET_OUTDIR / "files" / "feature_mapq_depth.tsv")
    if config["plotting"]["activate"]:
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_chrom_distribution.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_global_distribution.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_by_chrom.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "depth_by_regions.png",sample=SAMPLES))
        final_output.extend(expand(OUTDIR / "plots" / "{sample}" / "mapq.png",sample=SAMPLES))
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