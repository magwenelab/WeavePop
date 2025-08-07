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
    script:
        "../scripts/depth_by_windows.py"


# =================================================================================================
#   Per sample | Call CNVs
# =================================================================================================


rule cnv_calling:
    input:
        unpack(cnv_calling_input),
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
    script:
        "../scripts/cnv_calling.xsh"