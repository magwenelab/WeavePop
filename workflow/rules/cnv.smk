# =================================================================================================
#   Per sample | Normalize depth and by region
# =================================================================================================

rule depth_by_regions:
    input:
        depth = rules.mosdepth_good.output.bed,
        global_mode = rules.depth_distribution.output.global_mode
    output:
        OUTDIR / "depth_quality" / "{sample}" / "depth_by_regions.tsv"
    conda:
        "../envs/samtools.yaml"
    params:
        smoothing_size = config["cnv"]["smoothing_size"]
    log:
        "logs/samples/depth_quality/depth_by_regions_{sample}.log"
    shell:
        "xonsh workflow/scripts/depth_by_regions.xsh -di {input.depth} -gi {input.global_mode} -do {output} -s {params.smoothing_size} &> {log}"


# =================================================================================================
#   Per lineage | Run RepeatModeler and RepeatMasker
# =================================================================================================

rule repeat_modeler:
    input:
        rules.links.output
    output:
        known = REFDIR / "{lineage}" / "repeats" / "{lineage}_known.fa",
        unknown = REFDIR / "{lineage}" / "repeats" / "{lineage}_unknown.fa"
    params:
        repdir = "repeats"
    threads:
        config["cnv"]["repeats"]["repeats_threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/references/repeats/repeatmodeler_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-modeler.sh {threads} {input} {params.repdir} &> {log}"

rule repeat_masker:
    input:
        database = config["cnv"]["repeats"]["repeats_database"],
        fasta = rules.links.output,
        known = rules.repeat_modeler.output.known,
        unknown = rules.repeat_modeler.output.unknown
    output:
        REFDIR / "{lineage}" / "repeats" / "{lineage}_repeats.bed"
    threads:
        config["cnv"]["repeats"]["repeats_threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/references/repeats/repeatmasker_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-masker.sh {threads} {input.database} {input.fasta} {input.known} {input.unknown} {output} &> {log}"

# =================================================================================================
#   Per sample | Intercept depth by regions with repeats and call CNVs
# =================================================================================================

def cnv_calling_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "depth": OUTDIR / "depth_quality" / s["sample"] / "depth_by_regions.tsv",
        "repeats": REFDIR / s["lineage"]  / "repeats" / (s["lineage"] + "_repeats.bed")
        }
rule cnv_calling:
    input:
        unpack(cnv_calling_input)
    output:
        OUTDIR / "cnv" / "{sample}" / "cnv_calls.tsv"
    conda:
        "../envs/samtools.yaml"
    params:
        region_size = config["depth_quality"]["mosdepth"]["window"],
        depth_threshold = config["cnv"]["depth_threshold"]
    log:
        "logs/samples/cnv/cnv_calling_{sample}.log"
    shell:
        "xonsh workflow/scripts/cnv_calling.xsh -di {input.depth} -ri {input.repeats} -co {output} -sp {wildcards.sample} -rp {params.region_size} -dp {params.depth_threshold} &> {log}"
