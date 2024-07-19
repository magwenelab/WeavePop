# =================================================================================================
#   Join all lineages | Create table with loci to add to plots
# =================================================================================================

rule loci:
    input:
        refs = expand(REFDIR / "{lineage}" / "{lineage}.gff.tsv", lineage=LINEAGES),
        loci = LOCI_FILE
    output:
        locitable = DATASET_OUTDIR / "files"/ "loci_to_plot.tsv"
    log: 
        "logs/dataset/files/loci.log"
    shell:
        "xonsh workflow/scripts/loci.xsh -g {input.loci} -o {output} {input.refs} &> {log}"

# =================================================================================================
#   Per sample | Plot depth distribution and depth by chromosome 
# =================================================================================================

rule depth_distribution_plots:
    input:
        rules.depth_distribution.output.distrib,
        rules.depth_distribution.output.global_mode,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "depth_chrom_distribution.png",
        OUTDIR / "plots" / "{sample}" / "depth_global_distribution.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/depth_distribution_{sample}.log"
    script:
        "../scripts/depth_distribution_plots.R"

rule depth_by_chrom_plots:
    input:
        rules.depth_by_chrom_raw.output,
        rules.depth_by_chrom_good.output,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "depth_by_chrom.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/depth_by_chrom_{sample}.log"
    script:
        "../scripts/depth_by_chrom_plots.R"

# =================================================================================================
#   Per sample | Plot depth and mapq by regions
# =================================================================================================

def depth_by_regions_plots_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "depth": OUTDIR / "mosdepth" / s["sample"]  / "depth_by_regions.tsv",
        "cnv": OUTDIR / "cnv" / s["sample"] / "cnv_calls.tsv",
        "repeats": REFDIR / s["group"]  / "repeats" / (s["group"] + "_repeats.bed")
    }
rule depth_by_regions_plots:
    input:
        unpack(depth_by_regions_plots_input),
        loci = rules.loci.output.locitable,
        chrom_names = CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "depth_by_regions.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/depth_by_regions_{sample}.log"
    script:
        "../scripts/depth_by_regions_plots.R"

def mapq_plot_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "mapq": OUTDIR / "samtools" / s["sample"] / "mapq_window.bed",
        "structure": OUTDIR / "cnv" / s["sample"] / "cnv_calls.tsv",
        "repeats": REFDIR / s["group"]  / "repeats" / (s["group"] + "_repeats.bed")
    }
rule mapq_plot:
    input:
        unpack(mapq_plot_input),
        chrom_names = CHROM_NAMES,
        loci = rules.loci.output.locitable
    output:
        OUTDIR / "plots" / "{sample}" / "mapq.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/mapq_plot_{sample}.log"
    script:
        "../scripts/mapq.R"
