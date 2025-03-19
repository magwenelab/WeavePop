# =================================================================================================
#   Load modules
# =================================================================================================

import pandas as pd
import os.path
import glob
from pathlib import Path
import sys
from snakemake.utils import validate
import subprocess

# =================================================================================================
#   Print welcome message
# =================================================================================================

print("                                   ", flush=True)
print(" _           _   _      _   _   _  ", flush=True)
print("|_ | | |\\ | |_  |_| |  |_| | | |_|", flush=True)
print("|  |_| | \\| |_| | | |_ |   |_| |  ", flush=True)
print("                                   ", flush=True)
print("                                   ", flush=True)

# =================================================================================================
#   Print commit hash of current version
# =================================================================================================

def get_head_hash():
    try:
        result = subprocess.run(['sh', '-c', "tail -n1 .git/logs/HEAD | cut -d' ' -f2"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while getting the latest commit hash: {e.stderr}")
        return None
head_hash = get_head_hash()

if not head_hash:
    result = subprocess.run(['sh', '-c', "cat .head_hash"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    head_hash = result.stdout.strip()

print("Commmit hash of current version:")
print(f"{head_hash}")
print("", flush=True)

# =================================================================================================
#  Print configuration
# =================================================================================================

print("Executed command:", flush=True)
print(" ".join(sys.argv))
print("", flush=True)

profile_path = None
if "--profile" in sys.argv:
    profile_dir = sys.argv[sys.argv.index("--profile") + 1]
    profile_path = os.path.join(profile_dir, "config.yaml")
    print("", flush=True)
    print(".........................Execution Profile File:", profile_path, ".........................", flush=True)
    print("", flush=True)
    if profile_path and os.path.exists(profile_path):
        with open(profile_path, 'r') as file:
            print(file.read(), flush=True)

config_path = None
if "--configfile" in sys.argv:
    config_path = sys.argv[sys.argv.index("--configfile") + 1]
elif "--profile" in sys.argv:
    profile_dir = sys.argv[sys.argv.index("--profile") + 1]
    profile_path = os.path.join(profile_dir, "config.yaml")
    with open(profile_path, 'r') as file:
        for line in file:
            if line.startswith("configfile"):
                config_path = line.split(":")[1].strip()
                break
            else:
                config_path = "config/config.yaml"            
else:
    config_path = "config/config.yaml"

print("", flush=True)
print("............................Configuration File:", config_path, "............................", flush=True)
print("", flush=True)
if config_path and os.path.exists(config_path):
    with open(config_path, 'r') as file:
        print(file.read(), flush=True)

print("", flush=True)
print(".......................................Working directory.........................................", flush=True)
print("", flush=True)
print(os.getcwd(), flush=True)
print("", flush=True)
print("........................................Output directory.........................................", flush=True) 
print("", flush=True)
print(os.path.join(os.getcwd(), OUTPUT ), flush=True)
print("", flush=True)
print(".........................................Logs directory..........................................", flush=True) 
print("", flush=True)
print(os.path.join(os.getcwd(), LOGS ), flush=True)
print("", flush=True)

# =================================================================================================
#   Validate input files and get metadata table
# =================================================================================================

print("..............................................Input..............................................", flush=True)
print("                                   ", flush=True)

# --Define input data variables--------------------------------------------------------------------

if config["references_directory"].startswith("/"):
    REF_DATA = Path(config["references_directory"])
else:
    REF_DATA = Path(os.path.join(config["project_directory"], config["references_directory"]))

if config["fastqs_directory"].startswith("/"):
    FQ_DATA = Path(config["fastqs_directory"])
else:
    FQ_DATA = Path(os.path.join(config["project_directory"], config["fastqs_directory"]))

FQ1 = config["fastq_suffix1"]
FQ2 = config["fastq_suffix2"]

# --Validate original metadata table---------------------------------------------------------------

SAMPLE_ORIGINAL_FILE = Path(os.path.join(config["project_directory"], config["metadata"]))

if not SAMPLE_ORIGINAL_FILE:
    print("Metadata file not provided.", flush=True)
    print("Exiting...", flush=True)
    exit(1)
if os.path.exists(SAMPLE_ORIGINAL_FILE):
    print(f"Validating metadata file {SAMPLE_ORIGINAL_FILE}...", flush=True)
    SAMPLE_ORIGINAL = pd.read_csv(SAMPLE_ORIGINAL_FILE, sep=",", header=0)
    validate(SAMPLE_ORIGINAL, schema="../schemas/metadata.schema.yaml")
    print(f"    Number of samples in the metadata file: {SAMPLE_ORIGINAL.shape[0]}", flush=True)
    if SAMPLE_ORIGINAL.shape[0] == 0:
        print("Metadata file is empty.", flush=True)
        print("Exiting...", flush=True)
        exit(1)
else:
    print(f"Metadata file {SAMPLE_ORIGINAL_FILE} not found.", flush=True)
    print("Exiting...", flush=True)
    exit(1)

# --Exclude samples--------------------------------------------------------------------------------

if config["samples_to_exclude"]:
    EXCLUDE_FILE = Path(os.path.join(config["project_directory"], config["samples_to_exclude"]))    
    if os.path.exists(EXCLUDE_FILE):
        print(f"Excluding samples in the file {EXCLUDE_FILE} from the analysis...", flush=True)
        EXCLUDE_SAMPLES = set(list(pd.read_csv(EXCLUDE_FILE, header=None, names=["sample"])["sample"]))
        UNFILT_SAMPLE_TABLE = SAMPLE_ORIGINAL[~SAMPLE_ORIGINAL["sample"].isin(EXCLUDE_SAMPLES)]
        print(f"    Number of samples to analyze: {UNFILT_SAMPLE_TABLE.shape[0]}", flush=True)
    else:
        print(f"File {EXCLUDE_FILE} not found.", flush=True)
        print("Exiting...", flush=True)
        exit(1)
else:
    UNFILT_SAMPLE_TABLE = SAMPLE_ORIGINAL
    print(f"    Number of samples to analyze: {UNFILT_SAMPLE_TABLE.shape[0]}", flush=True)

# --Validate chromosome names file--------------------------------------------------------------------

print("                                  ", flush=True)
CHROM_NAMES = Path(os.path.join(config["project_directory"], config["chromosomes"]))

if not CHROM_NAMES:
    print("Chromosome names file not provided.", flush=True)
    print("Exiting...", flush=True)
    exit(1)
if os.path.exists(CHROM_NAMES):
    print(f"Validating chromosome names file {CHROM_NAMES}...", flush=True)
    CHROM_NAMES_TABLE = pd.read_csv(CHROM_NAMES, header=0, dtype={"chromosome": "string"})
    validate(CHROM_NAMES_TABLE, schema="../schemas/chromosomes.schema.yaml")
else:
    print(f"Chromosome names file {CHROM_NAMES} not found.", flush=True)
    print("Exiting...", flush=True)
    exit(1)

try:
    if all(item in CHROM_NAMES_TABLE["lineage"].unique() for item in UNFILT_SAMPLE_TABLE["lineage"].unique()):
        print("    All lineages are in the chromosome names file.", flush=True)
    else:
        print("Not all lineages from the metadata table are present in the chromosomes file.", flush=True)
        print("Exiting...", flush=True)
        exit(1)

    for lineage in UNFILT_SAMPLE_TABLE["lineage"].unique():
        print(f"Checking the chromosome names of lineage {lineage}...", flush=True)
        accessions = CHROM_NAMES_TABLE[CHROM_NAMES_TABLE["lineage"] == lineage]["accession"].tolist()
        ref_file = REF_DATA / f"{lineage}.fasta"
        if ref_file.exists():
            with open(ref_file) as f:
                seq_ids = [
                    line.strip().split()[0][1:] for line in f if line.startswith(">")
                ]

                if not all([acc in seq_ids for acc in accessions]):
                    raise ValueError(
                        f"Not all chromosomes in provided file are present in {ref_file}."
                    )
                else:
                    print(f"    All chromosomes in provided file are present in {ref_file}.", flush=True)
                                    
                if not all([seq_id in accessions for seq_id in seq_ids]):
                    raise ValueError(
                        f"Not all chromosomes in {ref_file} are present in the provided file."
                    )
                else:
                    print(f"    All chromosomes in {ref_file} are present in the provided file.", flush=True)

        else:
            raise ValueError(f"Reference {ref_file} not found.")

except Exception as e:
    print("Error in input files:", flush=True)
    print(e)
    print("Exiting...", flush=True)
    exit(1)

# --Validate loci file----------------------------------------------------------------------------

if config["plotting"]["loci"]:
    LOCI_FILE = Path(os.path.join(config["project_directory"], config["plotting"]["loci"]))
    if os.path.exists(LOCI_FILE):
        print("                                  ", flush=True)
        print("Validating provided file of loci to plot...", flush=True)
        LOCI_TABLE = pd.read_csv(LOCI_FILE, sep=",", header=0)
        validate(LOCI_TABLE, schema="../schemas/loci.schema.yaml")
        print(f"    Number of loci to plot: {LOCI_TABLE.shape[0]}", flush=True)
    else:
        print(f"File {LOCI_FILE} not found.", flush=True)
        print("Exiting...", flush=True)
        exit(1)
else:
    LOCI_FILE = SAMPLES_DIR / "loci_empty.txt"
    with open(LOCI_FILE, "w") as f:
        f.write("")
print("", flush=True)

# --Validate repeats file--------------------------------------------------------------------------

if config["cnv"]["activate"] or config["plotting"]["activate"] or config["database"]["activate"]:
    if config["cnv"]["repeats"]["repeats_database"]:
        if config["cnv"]["repeats"]["repeats_database"].startswith("/"):
            REPEATS_FILE = Path(config["cnv"]["repeats"]["repeats_database"])
            REPEATS_FILE = os.path.relpath(REPEATS_FILE, Path(os.getcwd()))
        else:
            REPEATS_FILE = Path(os.path.join(config["project_directory"], config["cnv"]["repeats"]["repeats_database"]))
        if not os.path.exists(REPEATS_FILE):
            print(f"Database of repetitive sequences file {REPEATS_FILE} not found.", flush=True)
            print("Exiting...", flush=True)
            exit(1)
        else:
            print(f"Database of repetitive sequences file {REPEATS_FILE} found.", flush=True)

    else:
        print("Database of repetitive sequences not provided.", flush=True)
        REPEATS_FILE = Path(os.path.join(config["project_directory"], "config/fake_repeats.fasta"))
        REPEATS_FILE.parent.mkdir(parents=True, exist_ok=True)
        with open(REPEATS_FILE, "w") as f:
            f.write(">fake\naaaaaaaaaaaaaa\n")
        print(f"WARNING: Using a fake database file {REPEATS_FILE}.", flush=True)
        print("The identification of repeats will not be accurate.", flush=True)
       


print("", flush=True)
print(".........................................Starting Workflow.........................................", flush=True)
print("", flush=True)

# =================================================================================================
#   Define directories and variables
# =================================================================================================

UNFILT_SAMPLES = list(set(UNFILT_SAMPLE_TABLE["sample"]))

# --Define directories----------------------------------------------------------------------------

SAMPLES_DIR = OUTPUT / SAMPLES_DIR_NAME
DATASET_DIR = OUTPUT / DATASET_DIR_NAME
REFS_DIR = OUTPUT / REFS_DIR_NAME
INTDIR = OUTPUT / INTDIR_NAME
INT_SAMPLES_DIR = INTDIR / SAMPLES_DIR_NAME
INT_DATASET_DIR = INTDIR / DATASET_DIR_NAME
INT_REFS_DIR = INTDIR / REFS_DIR_NAME
TEMPDIR = str(INTDIR / TEMPDIR_NAME)

# =================================================================================================
#   Variables for the module Reference annotation
# =================================================================================================

if config["annotate_references"]["activate"]:
    MAIN_FASTA = REF_DATA / config["annotate_references"]["fasta"]
    MAIN_GFF = REF_DATA / config["annotate_references"]["gff"]

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
        "mapq": INT_SAMPLES_DIR / "depth_quality" / s["sample"] / "mapq_by_window.bed",
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
    if config["depth_quality_features"]["activate"]:
        final_output.append(DATASET_DIR / "depth_quality" / "mapq_depth_by_feature.tsv")
    if config["snpeff"]["activate"]:
        final_output.append(DATASET_DIR / "snpeff" / "effects.tsv")
    if config["cnv"]["activate"]:
        final_output.append(DATASET_DIR / "cnv" / "cnv_calls.tsv")
    if config["plotting"]["activate"]:
        final_output.append(DATASET_DIR / "plots" / "dataset_depth_by_chrom.png")
        final_output.append(DATASET_DIR / "plots" / "dataset_summary.png")
    if config["database"]["activate"]:
        final_output.append(expand(DATASET_DIR / "database.db"))
    return final_output
