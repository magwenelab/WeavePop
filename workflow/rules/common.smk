import pandas as pd
import os.path
import glob

#### Defining global variables ####
samplefile=(pd.read_csv(config["sample_file"], sep=","))
SAMPLES=list(set(samplefile["sample"]))
LINEAGES=list(set(samplefile["group"]))
REF_DATA = config["references"]["directory"]
# REF_FASTAS = set(glob.glob(f"{REF_DIR}/*.fasta") + glob.glob(f"{REF_DIR}/*.fa") + glob.glob(f"{REF_DIR}/*.fna") + glob.glob(f"{REF_DIR}/*.fas"))
# REF_NAMES = [os.path.splitext(os.path.basename(f))[0] for f in REF_FASTAS]
OUTDIR= "results/samples"
DATASET_OUTDIR = "results/dataset"
REFDIR = "results/references"

FEATURE_FILE = "config/features.txt"

#### Defining variables for the reference annotation module(references.smk) ####
if config["annotate_references"]["activate"]:
    MAIN_DIR = config["annotate_references"]["directory"]
    MAIN_FASTA = os.path.join(MAIN_DIR, config["annotate_references"]["fasta"])
    MAIN_GFF = os.path.join(MAIN_DIR, config["annotate_references"]["gff"])
    MAIN_NAME, _ = os.path.splitext(os.path.basename(MAIN_GFF))
#### Defining which final output files are being requested ####
def get_final_output():
    # final_output = expand("analysis/{sample}/snps.bam",sample=SAMPLES)
    # final_output.extend(expand("analysis/{sample}/lifted.gff_polished", sample=SAMPLES))

    if config["annotate_references"]["activate"]:
        final_output = (expand(os.path.join(REFDIR, "{lineage}.gff"),lineage=LINEAGES))
        final_output.append(str(os.path.join(REFDIR, MAIN_NAME + ".tsv")))
    return final_output

#### Creating softlinks to have the reference genomes in the REFDIR ####
rule links:
    input:
        os.path.join(REF_DATA, "{lineage}.fasta")
    output:
        os.path.join(REFDIR, "{lineage}.fasta")
    shell:
        "ln -s -r {input} {output}"
