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

OUTPUT = Path(config["joint_output_directory"])
SAMPLES_DIR = OUTPUT / SAMPLES_DIR_NAME
DATASET_DIR = OUTPUT / DATASET_DIR_NAME
REFS_DIR = OUTPUT / REFS_DIR_NAME
INTDIR = OUTPUT / INTDIR_NAME
INT_SAMPLES_DIR = INTDIR / SAMPLES_DIR_NAME
INT_DATASET_DIR = INTDIR / DATASET_DIR_NAME
INT_REFS_DIR = INTDIR / REFS_DIR_NAME

INPUT_PATHS = config["datasets_paths"].split(",")
LIST_PATHS = [Path(dir) for dir in INPUT_PATHS]

# =================================================================================================
#   Define final output
# =================================================================================================

def get_final_output():
    final_output = DATASET_DIR / "database.db"
    return final_output

# =================================================================================================
#   Join the metadata, chromosomes and loci files from the different datasets
# =================================================================================================

rule join_results:
    input:
        metadata_files = expand(os.path.join("{dir}", INTDIR_NAME, "metadata.csv"), dir = LIST_PATHS),
        chromosomes_files = expand(os.path.join("{dir}", INTDIR_NAME, "chromosomes.csv"), dir = LIST_PATHS),
        loci_files = expand(os.path.join("{dir}", INTDIR_NAME, "loci.csv"), dir = LIST_PATHS)
    output:
        metadata = INTDIR / "metadata.csv",
        chromosomes = INTDIR / "chromosomes.csv",
        loci = INTDIR / "loci.csv"
    params:
        dataset_names = config["datasets_names"].split(",")
    run:
        dataframes = []
        for file_path, string in zip(input.metadata_files, params.dataset_names):
            df = pd.read_csv(file_path, sep=",", header=0)
            df["dataset"] = string
            dataframes.append(df)

        metadata= pd.concat(dataframes, ignore_index=True)
        metadata = metadata.drop_duplicates()
        chromosomes = pd.concat([pd.read_csv(f, sep=",", header=0) for f in input.chromosomes_files])
        chromosomes = chromosomes.drop_duplicates()
        loci = pd.concat([pd.read_csv(f, sep=",", header=0) for f in input.loci_files])
        loci = loci.drop_duplicates()
        metadata.to_csv(str(output.metadata), sep=",", index=False)
        chromosomes.to_csv(str(output.chromosomes), sep=",", index=False)
        loci.to_csv(str(output.loci), sep=",", index=False)

# Make the list of lineages from the joined metadata
checkpoint make_lineages_dirs:
    input:
        metadata = rules.join_results.output.metadata
    output:
        directory(INT_REFS_DIR / "lineage_names")
    run:
        metadata = pd.read_csv(input.metadata, sep=",", header=0)
        lineages = list(set(metadata["lineage"]))
        for lineage in lineages:
            path = Path( INT_REFS_DIR / "lineage_names") / f"{lineage}.lineage"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()

def listing_lineages(wildcards):
    checkpoint_output = checkpoints.make_lineages_dirs.get(**wildcards).output[0]
    return expand("{i}",
                i=glob_wildcards(os.path.join(checkpoint_output, "{i}.lineage" )).i)
LINEAGES = listing_lineages


# =================================================================================================
#   Join tabular results
# =================================================================================================
def input_joining(wildcards):
    paths_cnv = []
    paths_mapq_depth = []
    paths_cds = []
    paths_prots = []
    for dir in LIST_PATHS:
        metadata = os.path.join(dir, INTDIR_NAME, "metadata.csv")
        samps_dir_df = pd.read_csv(metadata, header=0)
        samps = samps_dir_df["sample"]
        for samp in samps:
            cnv = os.path.join(dir, SAMPLES_DIR_NAME, "cnv", samp, "cnv_calls.tsv")
            mapq_depth = os.path.join(dir, SAMPLES_DIR_NAME, "depth_quality", samp, "feature_mapq_depth.tsv")
            cds = os.path.join(dir, INTDIR_NAME, SAMPLES_DIR_NAME, "annotation", samp, "cds.csv")
            prots = os.path.join(dir, INTDIR_NAME, SAMPLES_DIR_NAME, "annotation", samp, "proteins.csv")
            paths_cnv.append(cnv)
            paths_mapq_depth.append(mapq_depth)
            paths_cds.append(cds)
            paths_prots.append(prots)
    return {
        "cnv": paths_cnv,
        "mapq_depth": paths_mapq_depth,
        "cds": paths_cds,
        "prots": paths_prots
    }

def get_input_cnv(wildcards):
    return input_joining(wildcards)["cnv"]

def get_input_mapq_depth(wildcards):
    return input_joining(wildcards)["mapq_depth"]

def get_input_cds(wildcards):
    return input_joining(wildcards)["cds"]

def get_input_prots(wildcards):
    return input_joining(wildcards)["prots"]

rule join_sequences:
    input:
        cds=get_input_cds,
        prots=get_input_prots
    output:
        sequences = INT_DATASET_DIR / "sequences.csv"
    run:
        cds = pd.concat([pd.read_csv(f, sep="\t") for f in input.cds])
        prots = pd.concat([pd.read_csv(f, sep="\t") for f in input.prots])
        sequences = pd.concat([cds, prots])
        sequences.to_csv(output.sequences, sep="\t", index=False)

rule join_cnv:
    input:
        cnv=get_input_cnv,
    output:
        cnv = DATASET_DIR / "cnv" / "cnv_calls.tsv"
    run:
        cnv = pd.concat([pd.read_csv(f, sep="\t") for f in input.cnv])
        cnv.to_csv(output.cnv, sep="\t", index=False)

rule join_mapq_depth:
    input:
        mapq_depth=get_input_mapq_depth,
    output:
        mapq_depth = DATASET_DIR / "depth_quality" / "feature_mapq_depth.tsv",
    run: 
        mapq_depth = pd.concat([pd.read_csv(f, sep="\t") for f in input.mapq_depth])
        mapq_depth.to_csv(output.mapq_depth, sep="\t", index=False)

def input_gffs(wildcards):
    paths = []
    for dir in LIST_PATHS:
        metadata = os.path.join(dir, INTDIR_NAME, "metadata.csv")
        lineages_dir_df = pd.read_csv(metadata, header=0)
        lineages = set(lineages_dir_df["lineage"])
        for lineage in lineages:
            gff = os.path.join(dir, INTDIR_NAME, REFS_DIR_NAME, lineage, f"{lineage}.gff.tsv")
            paths.append(gff)
    return paths

rule join_gffs:
    input:
        input_gffs
    output:
        INT_REFS_DIR / "all_lineages.gff.tsv"
    log: 
        "logs/references/join_gffs.log"
    shell:
        "python workflow/scripts/join_gffs.py -o {output} {input} &> {log}"

# =================================================================================================
#   Copy SnpEff lineage databases and config file
# =================================================================================================

def input_copy_speff_data(wildcards):
    paths_snpeff_Data = {}
    for dir in LIST_PATHS:
        metadata = os.path.join(dir, INTDIR_NAME, "metadata.csv")
        lins_dir_df = pd.read_csv(metadata, header=0)
        lins = lins_dir_df["lineage"]
        for lin in lins:
            data_dict = {lin: os.path.join(dir, INTDIR_NAME, REFS_DIR_NAME, "snpeff_data", config["species_name"] + f"_{lin}")}
            paths_snpeff_Data.update(data_dict)
    list_paths = list(paths_snpeff_Data.values())
    return list_paths

rule copy_snpeff_data:
    input:
        input_copy_speff_data
    output:
        INT_REFS_DIR / "snpeff_data" / "copy.done"
    params:
        dir = INT_REFS_DIR / "snpeff_data"
    shell:
        """
        ln -srf {input} {params.dir} && touch {output}
        """

rule copy_snpeff_config:
    input:
        data = rules.copy_snpeff_data.output,
        config = expand(os.path.join("{dir}", INTDIR_NAME, REFS_DIR_NAME, "snpeff_data", "snpEff.config"), dir = LIST_PATHS)
    output:
        INT_REFS_DIR / "snpeff_data" / "snpEff.config"
    shell:
        "cat {input.config} > {output}"

# =================================================================================================
#   Intersect SNPs, run SnpEff and extract annotations
# =================================================================================================

def intersect_vcfs_input(wildcards):
    paths = []
    for dir in LIST_PATHS:
        metadata = os.path.join(dir, INTDIR_NAME, "metadata.csv")
        samps_dir_df = pd.read_csv(metadata, header=0)
        samps_dir_df = samps_dir_df.loc[samps_dir_df["lineage"] == wildcards.lineage]
        samps = samps_dir_df["sample"]
        for samp in samps:
            vcf = os.path.join(dir, SAMPLES_DIR_NAME, "snippy", samp , "snps.vcf.gz")
            paths.append(vcf)
    return {
        "vcfs" : paths
    }

rule intersect_vcfs:
    input:
        unpack(intersect_vcfs_input)
    output:
        vcf = INT_DATASET_DIR / "snps" / "{lineage}_intersection.vcf",
        tsv = INT_DATASET_DIR / "snps" / "{lineage}_presence.tsv"
    params:
        tmp_dir = "tmp_{lineage}"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/intersect_vcfs_{lineage}.log"
    shell:
        "xonsh workflow/scripts/intersect_vcfs.xsh "
        "-v {output.vcf} "
        "-p {output.tsv} "
        "-l {wildcards.lineage} "
        "-t {params.tmp_dir} "
        "{input.vcfs} "
        "&> {log}"

rule snpeff:
    input:
        vcf = rules.intersect_vcfs.output.vcf,
        config = rules.copy_snpeff_config.output
    output:
        vcf = INT_DATASET_DIR / "snps" / "{lineage}_snpeff.vcf",
        html = INT_DATASET_DIR / "snps" / "{lineage}_snpeff.html"
    params:
        dir = os.getcwd() / INT_REFS_DIR / "snpeff_data",
        name = config["species_name"] + "_{lineage}"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/snpeff_{lineage}.log"
    shell:
        "snpEff ann -v -classic "
        "-dataDir {params.dir} "
        "-config {input.config} "
        "-s {output.html} "
        "{params.name} "
        "{input.vcf} "
        "1> {output.vcf} 2> {log}"

rule extract_vcf_annotation:
    input:
        vcf = rules.snpeff.output.vcf
    output:
        effects = INT_DATASET_DIR / "snps" / "{lineage}_effects.tsv",
        variants = INT_DATASET_DIR / "snps" / "{lineage}_variants.tsv",
        lofs = INT_DATASET_DIR / "snps" / "{lineage}_lofs.tsv",
        nmds = INT_DATASET_DIR / "snps" / "{lineage}_nmds.tsv"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/extract_vcf_annotation_{lineage}.log"
    shell:
        "xonsh workflow/scripts/extract_vcf_annotation.xsh "
        "-i {input.vcf} "
        "-e {output.effects} "
        "-v {output.variants} "
        "-f {output.lofs} "
        "-n {output.nmds} "
        "-l {wildcards.lineage} "
        "&> {log}"
    
# =================================================================================================
#   Create final database
# =================================================================================================
rule join_variant_annotation:
    input:
        effects = expand(INT_DATASET_DIR / "snps" / "{lineage}_effects.tsv", lineage=LINEAGES),
        variants = expand(INT_DATASET_DIR / "snps" / "{lineage}_variants.tsv", lineage=LINEAGES),
        lofs = expand(INT_DATASET_DIR / "snps" / "{lineage}_lofs.tsv", lineage=LINEAGES),
        nmds = expand(INT_DATASET_DIR / "snps" / "{lineage}_nmds.tsv", lineage=LINEAGES),
        presence = expand(INT_DATASET_DIR / "snps" / "{lineage}_presence.tsv", lineage=LINEAGES)
    output:
        effects = DATASET_DIR / "snps" / "effects.tsv",
        variants = DATASET_DIR / "snps" / "variants.tsv",
        lofs = DATASET_DIR / "snps" / "lofs.tsv",
        nmds = DATASET_DIR / "snps" / "nmds.tsv",
        presence = DATASET_DIR / "snps" / "presence.tsv"
    run:
        effects = pd.concat([pd.read_csv(f, sep="\t") for f in input.effects])
        variants = pd.concat([pd.read_csv(f, sep="\t") for f in input.variants])
        lofs = pd.concat([pd.read_csv(f, sep="\t") for f in input.lofs])
        nmds = pd.concat([pd.read_csv(f, sep="\t") for f in input.nmds])
        presence = pd.concat([pd.read_csv(f, sep="\t") for f in input.presence])
        effects.to_csv(output.effects, sep = "\t", index=False)
        variants.to_csv(output.variants, sep="\t", index=False)
        lofs.to_csv(output.lofs, sep="\t", index=False)
        nmds.to_csv(output.nmds, sep="\t", index=False)
        presence.to_csv(output.presence, sep="\t", index=False)

rule complete_db:
    input:
        metadata = INTDIR / "metadata.csv",
        chrom_names = INTDIR / "chromosomes.csv",
        cnv = rules.join_cnv.output.cnv,
        md = rules.join_mapq_depth.output.mapq_depth,
        gffs = rules.join_gffs.output,
        effects = rules.join_variant_annotation.output.effects,
        variants = rules.join_variant_annotation.output.variants,
        presence = rules.join_variant_annotation.output.presence,
        lofs = rules.join_variant_annotation.output.lofs,
        nmds = rules.join_variant_annotation.output.nmds,
        seqs = rules.join_sequences.output.sequences
    output:
        DATASET_DIR / "database.db"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/complete_db.log"
    shell:
        "xonsh workflow/scripts/build_database.xsh "
        "-m {input.metadata} "
        "-ch {input.chrom_names} "
        "-cnv {input.cnv} "
        "-md {input.md} "
        "-g {input.gffs} "
        "-e {input.effects} "
        "-v {input.variants} "
        "-p {input.presence} "
        "-l {input.lofs} "
        "-n {input.nmds} "
        "-s {input.seqs} "
        "-o {output} &> {log}"

