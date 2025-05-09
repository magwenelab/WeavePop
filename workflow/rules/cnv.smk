# =================================================================================================
#   Per sample | Run Mosdepth to get depth per window of good quality reads
# =================================================================================================


rule mosdepth:
    input:
        bam=INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "snps_good.bam",
        bai=INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "snps_good.bam.bai",
    output:
        bed=INT_SAMPLES_DIR / "mosdepth" / "{sample}" / "coverage_good.regions.bed.gz",
    params:
        window=config["depth_quality"]["mosdepth"]["window"],
        extra=config["depth_quality"]["mosdepth"]["extra"],
        min_mapq=config["depth_quality"]["flag_quality"]["min_mapq"],
        outdir=INT_SAMPLES_DIR / "mosdepth",
    log:
        LOGS / "samples" / "depth_quality" / "mosdepth_{sample}.log",
    threads: config["depth_quality"]["mosdepth"]["threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/depth.yaml"
    shell:
        "mosdepth "
        "-n "
        "--by {params.window} "
        "--mapq {params.min_mapq} "
        "-t {threads} "
        "{params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage_good "
        "{input.bam} "
        "&> {log}"


# =================================================================================================
#   Per sample | Normalize and smooth depth by windows
# =================================================================================================


rule depth_by_windows:
    input:
        depth=rules.mosdepth.output.bed,
    output:
        windows=INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_windows.tsv",
        chroms=SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_chrom.tsv"
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
        "-co {output.chroms} "
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
