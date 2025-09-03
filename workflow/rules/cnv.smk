# =================================================================================================
#   Per sample | Call CNVs
# =================================================================================================


rule cnv_calling:
    input:
        unpack(cnv_calling_input),
    output:
        windows=INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_windows.tsv",
        cnvs=SAMPLES_DIR / "cnv" / "{sample}" / "cnv_calls.tsv",
        chrom_cnvs=SAMPLES_DIR / "cnv" / "{sample}" / "cnv_chromosomes.tsv",
    params:
        smoothing_size=config["cnv"]["smoothing_size"],
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
