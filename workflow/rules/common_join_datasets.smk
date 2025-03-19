# =================================================================================================
#   Load modules
# =================================================================================================

import pandas as pd
import os.path
import glob
from pathlib import Path
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
print(os.path.join(os.getcwd(), OUTPUT), flush=True)
print("", flush=True)
print(".........................................Logs directory..........................................", flush=True) 
print("", flush=True)
print(os.path.join(os.getcwd(), LOGS), flush=True)
print("", flush=True)
# =================================================================================================
#   Define global variables
# =================================================================================================

SAMPLES_DIR = OUTPUT / SAMPLES_DIR_NAME
DATASET_DIR = OUTPUT / DATASET_DIR_NAME
REFS_DIR = OUTPUT / REFS_DIR_NAME
INTDIR = OUTPUT / INTDIR_NAME
INT_SAMPLES_DIR = INTDIR / SAMPLES_DIR_NAME
INT_DATASET_DIR = INTDIR / DATASET_DIR_NAME
INT_REFS_DIR = INTDIR / REFS_DIR_NAME
TEMPDIR = str(INTDIR / TEMPDIR_NAME)

INPUT_PATHS = config["datasets_paths"].split(",")
LIST_PATHS = [Path(dir) for dir in INPUT_PATHS]

# =================================================================================================
#   Validate input files
# =================================================================================================

print("..............................................Input..............................................", flush=True)
print("                                   ", flush=True)

if OUTPUT in LIST_PATHS:
    print("", flush=True)
    print("Output directory is the same as an input directory. Try adding a run_id.", flush=True)
    print("Exiting...", flush=True)
    exit(1)

for dir in LIST_PATHS:
    if os.path.exists(dir):
        print(f"Dataset directory: {dir}")
    else:
        print("", flush=True)
        print(f"Directory {dir} not found.", flush=True)
        print("Exiting...", flush=True)
        exit(1)
    if os.path.exists(os.path.join(dir, DATASET_DIR_NAME, "metadata.csv")):
        print(f"Metadata file: {os.path.join(dir, DATASET_DIR_NAME, 'metadata.csv')}")
    else:
        print("", flush=True)
        print(f"File {os.path.join(dir, DATASET_DIR_NAME, 'metadata.csv')} not found", flush=True)
        print("Exiting...", flush=True)
        exit(1)
    if os.path.exists(os.path.join(dir, DATASET_DIR_NAME, "chromosomes.csv")):
        print(f"Chromosomes file: {os.path.join(dir, DATASET_DIR_NAME, 'chromosomes.csv')}")
    else:
        print("", flush=True)
        print(f"File {os.path.join(dir, DATASET_DIR_NAME, 'chromosomes.csv')} not found", flush=True)
        print("Exiting...", flush=True)
        exit(1)
print("", flush=True)
print(".........................................Starting Workflow.........................................", flush=True)
print("", flush=True)

# =================================================================================================
#   Input functions for rules
# =================================================================================================

def input_joining(wildcards):
    paths_cnv = []
    paths_mapq_depth = []
    paths_cds = []
    paths_prots = []
    for dir in LIST_PATHS:
        metadata = os.path.join(dir, DATASET_DIR_NAME, "metadata.csv")
        samps_dir_df = pd.read_csv(metadata, header=0)
        samps = samps_dir_df["sample"]
        for samp in samps:
            cnv = os.path.join(dir, SAMPLES_DIR_NAME, "cnv", samp, "cnv_calls.tsv")
            mapq_depth = os.path.join(
                dir, SAMPLES_DIR_NAME, "depth_quality", samp, "mapq_depth_by_feature.tsv"
            )
            cds = os.path.join(dir, INTDIR_NAME, SAMPLES_DIR_NAME, "annotation", samp, "cds.csv")
            prots = os.path.join(
                dir, INTDIR_NAME, SAMPLES_DIR_NAME, "annotation", samp, "proteins.csv"
            )
            paths_cnv.append(cnv)
            paths_mapq_depth.append(mapq_depth)
            paths_cds.append(cds)
            paths_prots.append(prots)
    return {
        "cnv": paths_cnv,
        "mapq_depth": paths_mapq_depth,
        "cds": paths_cds,
        "prots": paths_prots,
    }


def input_join_cnv(wildcards):
    return input_joining(wildcards)["cnv"]


def input_join_mapq_depth(wildcards):
    return input_joining(wildcards)["mapq_depth"]


def input_join_cds(wildcards):
    return input_joining(wildcards)["cds"]


def input_join_prots(wildcards):
    return input_joining(wildcards)["prots"]


def input_join_ref_annotations(wildcards):
    paths = []
    for dir in LIST_PATHS:
        metadata = os.path.join(dir, DATASET_DIR_NAME, "metadata.csv")
        lineages_dir_df = pd.read_csv(metadata, header=0)
        lineages = set(lineages_dir_df["lineage"])
        for lineage in lineages:
            gff = os.path.join(dir, REFS_DIR_NAME, lineage, f"{lineage}.gff.tsv")
            paths.append(gff)
    return paths


def input_copy_speff_data(wildcards):
    paths_snpeff_Data = {}
    for dir in LIST_PATHS:
        metadata = os.path.join(dir, DATASET_DIR_NAME, "metadata.csv")
        lins_dir_df = pd.read_csv(metadata, header=0)
        lins = lins_dir_df["lineage"]
        for lin in lins:
            data_dict = {
                lin: os.path.join(
                    dir,
                    INTDIR_NAME,
                    REFS_DIR_NAME,
                    "snpeff_data",
                    f"Species_name_{lin}",
                )
            }
            paths_snpeff_Data.update(data_dict)
    list_paths = list(paths_snpeff_Data.values())
    return list_paths


def input_intersect_vcfs(wildcards):
    paths = []
    for dir in LIST_PATHS:
        metadata = os.path.join(dir, DATASET_DIR_NAME, "metadata.csv")
        samps_dir_df = pd.read_csv(metadata, header=0)
        samps_dir_df = samps_dir_df.loc[samps_dir_df["lineage"] == wildcards.lineage]
        samps = samps_dir_df["sample"]
        for samp in samps:
            vcf = os.path.join(dir, SAMPLES_DIR_NAME, "snippy", samp, "snps.vcf.gz")
            paths.append(vcf)
    return {"vcfs": paths}


def input_symlink_ref_gff(wildcards):
    paths = {}
    for dir in LIST_PATHS:
        metadata = os.path.join(dir, DATASET_DIR_NAME, "metadata.csv")
        lineages_dir_df = pd.read_csv(metadata, header=0)
        lineages = set(lineages_dir_df["lineage"])
        for lineage in lineages:
            gff = os.path.join(dir, INTDIR_NAME, REFS_DIR_NAME, lineage, f"{lineage}_interg_introns.gff.tsv")
            paths[lineage] = gff
    lineage_path = paths[wildcards.lineage]
    return lineage_path

# =================================================================================================
#   Checkpoint functions
# =================================================================================================

def listing_lineages(wildcards):
    checkpoint_output = checkpoints.get_lineages.get(**wildcards).output[0]
    return expand("{i}", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.lineage")).i)


LINEAGES = listing_lineages

# =================================================================================================
#   Final output definition functions
# =================================================================================================

def get_final_output():
    final_output = DATASET_DIR / "database.db"
    return final_output
