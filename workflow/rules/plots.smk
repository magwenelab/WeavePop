# =================================================================================================
#   Join all lineages | Create table with loci to add to plots
# =================================================================================================

rule loci:
    input:
        refs = expand(REFDIR / "{lineage}" / "{lineage}.gff.tsv", lineage=LINEAGES),
        loci = rules.copy_config.output.l
    output:
        locitable = REFDIR / "loci_to_plot.tsv"
    log: 
        "logs/dataset/files/loci.log"
    shell:
        "xonsh workflow/scripts/loci.xsh -g {input.loci} -o {output} {input.refs} &> {log}"

# =================================================================================================
#   Per sample | Plot depth distribution and depth by chromosome 
# =================================================================================================

rule depth_distribution_plots:
    input:
        OUTDIR / "depth_quality" / "{sample}" / "depth_distribution.tsv",
        rules.copy_config.output.c
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
        OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_raw.tsv",
        OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv",
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
#   Per sample | Plot depth and mapq by windows
# =================================================================================================

def depth_by_windows_plots_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "depth": OUTDIR / "depth_quality" / s["sample"]  / "depth_by_windows.tsv",
        "cnv": OUTDIR / "cnv" / s["sample"] / "cnv_calls.tsv",
        "repeats": REFDIR / s["lineage"]  / "repeats" / (s["lineage"] + "_repeats.bed")
    }
rule depth_by_windows_plots:
    input:
        unpack(depth_by_windows_plots_input),
        loci = rules.loci.output.locitable,
        chrom_names = CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "depth_by_windows.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/depth_by_windows_{sample}.log"
    script:
        "../scripts/depth_by_windows_plots.R"

def mapq_plot_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "mapq": OUTDIR / "depth_quality" / s["sample"] / "mapq_window.bed",
        "cnv": OUTDIR / "cnv" / s["sample"] / "cnv_calls.tsv",
        "repeats": REFDIR / s["lineage"]  / "repeats" / (s["lineage"] + "_repeats.bed")
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
