# =================================================================================================
#   Per sample | Normalize and smooth depth by windows
# =================================================================================================


rule depth_by_windows:
    input:
        depth=rules.mosdepth.output.bed,
    output:
        windows=INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_windows.tsv",
    params:
        smoothing_size=config["cnv"]["smoothing_size"],
    log:
        LOGS / "samples" / "depth_quality" / "depth_by_windows_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/depth_by_windows.xsh "
        "-di {input.depth} "
        "-do {output.windows} "
        "-ss {params.smoothing_size} "
        "-sm {wildcards.sample} "
        "&> {log}"


# =================================================================================================
#   Per sample | Call CNVs
# =================================================================================================


rule cnv_calling:
    input:
        unpack(cnv_calling_input),
        chrom_length=rules.quality_filter.output.chromosomes,
    output:
        cnvs=SAMPLES_DIR / "cnv" / "{sample}" / "cnv_calls.tsv",
        chrom_cnvs=SAMPLES_DIR / "cnv" / "{sample}" / "cnv_chromosomes.tsv",
    params:
        window_size=config["depth_quality"]["mosdepth"]["window"],
        depth_threshold=config["cnv"]["depth_threshold"],
    log:
        LOGS / "samples" / "cnv" / "cnv_calling_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/cnv_calling.xsh "
        "-di {input.depth} "
        "-ri {input.repeats} "
        "-ai {input.annotation} "
        "-ci {input.chrom_length} "
        "-co {output.cnvs} "
        "-mo {output.chrom_cnvs} "
        "-sp {wildcards.sample} "
        "-wp {params.window_size} "
        "-dp {params.depth_threshold} "
        "-t {resources.tmpdir} &> {log}"
