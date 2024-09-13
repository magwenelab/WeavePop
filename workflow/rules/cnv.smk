# =================================================================================================
#   Per sample | Run Mosdepth to get depth per window of  good quality reads
# =================================================================================================

rule mosdepth_good:
    input:
        bam = SAMPLES_DIR / "snippy" / "{sample}" / "snps.bam",
        bai = SAMPLES_DIR / "snippy" / "{sample}" / "snps.bam.bai"
    output:
        bed = INT_SAMPLES_DIR / "mosdepth" / "{sample}" / "coverage_good.regions.bed.gz"
    params:
        window = config["depth_quality"]["mosdepth"]["window"],
        extra = config["depth_quality"]["mosdepth"]["extra"],
        min_mapq = config["depth_quality"]["mosdepth"]["min_mapq"],
        outdir = INT_SAMPLES_DIR / "mosdepth"
    log:
        "logs/samples/mosdepth/mosdepth_good_{sample}.log"
    threads:
       config["depth_quality"]["mosdepth"]["threads"] 
    resources:
        tmpdir = TEMPDIR
    conda: 
        "../envs/depth.yaml"
    shell:
        "mosdepth -n --by {params.window} --mapq {params.min_mapq} -t {threads} {params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage_good {input.bam} "
        "&> {log}"

# =================================================================================================
#   Per sample | Normalize depth and by windows
# =================================================================================================

rule depth_by_windows:
    input:
        depth = rules.mosdepth_good.output.bed,
        global_mode = SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv"
    output:
        INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_windows.tsv"
    params:
        smoothing_size = config["cnv"]["smoothing_size"]
    log:
        "logs/samples/depth_quality/depth_by_windows_{sample}.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/depth_by_windows.xsh -di {input.depth} -gi {input.global_mode} -do {output} -s {params.smoothing_size} &> {log}"

# =================================================================================================
#   Per lineage | Run RepeatModeler and RepeatMasker
# =================================================================================================

rule repeat_modeler:
    input:
        rules.links.output
    output:
        known = INT_REFS_DIR / "{lineage}" / "repeats" / "{lineage}_known.fa",
        unknown = INT_REFS_DIR / "{lineage}" / "repeats" / "{lineage}_unknown.fa"
    params:
        repdir = "repeats"
    log:
        "logs/references/repeats/repeatmodeler_{lineage}.log"
    threads:
        config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "bash workflow/scripts/repeat-modeler.sh {threads} {input} {params.repdir} &> {log}"

rule repeat_masker_1:
    input:
        database = config["cnv"]["repeats"]["repeats_database"],
        fasta = rules.links.output
    output: 
        INT_REFS_DIR / "{lineage}" / "repeats" / "01_simple" / "{lineage}.fasta.out"
    log:
        "logs/references/repeats/repeatmasker1_{lineage}.log"
    threads:
        config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        outdir=$(dirname {output})
        RepeatMasker -pa {threads} -lib {input.database} -a -e ncbi -dir $outdir -noint -xsmall {input.fasta} &> {log}
        """

rule repeat_masker_2:
    input:
        database = config["cnv"]["repeats"]["repeats_database"],
        fasta = rules.links.output
    output: 
        INT_REFS_DIR / "{lineage}" / "repeats" / "02_complex" / "{lineage}.fasta.out"
    log:
        "logs/references/repeats/repeatmasker2_{lineage}.log"
    threads:
        config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        outdir=$(dirname {output})
        RepeatMasker -pa {threads} -lib {input.database} -a -e ncbi -dir $outdir -nolow {input.fasta} &> {log}
        """

rule repeat_masker_3:
    input:
        known = rules.repeat_modeler.output.known,
        fasta = rules.links.output
    output: 
        INT_REFS_DIR / "{lineage}" / "repeats" / "03_known" / "{lineage}.fasta.out"
    log:
        "logs/references/repeats/repeatmasker3_{lineage}.log"
    threads:
        config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        outdir=$(dirname {output})
        RepeatMasker -pa {threads} -lib {input.known} -a -e ncbi -dir $outdir -nolow {input.fasta} &> {log}
        """

rule repeat_masker_4:
    input:
        unknown = rules.repeat_modeler.output.unknown,
        fasta = rules.links.output
    output: 
        INT_REFS_DIR / "{lineage}" / "repeats" / "04_unknown" / "{lineage}.fasta.out"
    log:
        "logs/references/repeats/repeatmasker4_{lineage}.log"
    threads:
        config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        outdir=$(dirname {output})
        RepeatMasker -pa {threads} -lib {input.unknown} -a -e ncbi -dir $outdir -nolow {input.fasta} &> {log}
        """

rule repeat_masker_bed:
    input:
        simple = rules.repeat_masker_1.output,
        complx = rules.repeat_masker_2.output,
        known = rules.repeat_masker_3.output,
        unknown = rules.repeat_masker_4.output
    output:
        simple = INT_REFS_DIR / "{lineage}" / "repeats" / "01_simple" / "{lineage}.bed",
        complx = INT_REFS_DIR / "{lineage}" / "repeats" / "02_complex" / "{lineage}.bed",
        known = INT_REFS_DIR / "{lineage}" / "repeats" / "03_known" / "{lineage}.bed",
        unknown = INT_REFS_DIR / "{lineage}" / "repeats" / "04_unknown" / "{lineage}.bed"
    log:
        "logs/references/repeats/repeatmasker_combine_{lineage}.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        tail -n +4 {input.simple} | awk '{{print $5"\t"($6-1)"\t"$7"\t"$11}}' 1> {output.simple} 2> {log}
        tail -n +4 {input.complx} | awk '{{print $5"\t"($6-1)"\t"$7"\t"$11}}' 1> {output.complx} 2>> {log}
        tail -n +4 {input.known} | awk '{{print $5"\t"($6-1)"\t"$7"\t"$11}}' 1> {output.known} 2>> {log}
        tail -n +4 {input.unknown} | awk '{{print $5"\t"($6-1)"\t"$7"\t"$11}}' 1> {output.unknown} 2>> {log}
        """
    
rule repeat_masker_combine:
    input:
        simple = rules.repeat_masker_bed.output.simple,
        complx = rules.repeat_masker_bed.output.complx,
        known = rules.repeat_masker_bed.output.known,
        unknown = rules.repeat_masker_bed.output.unknown
    output:
        REFS_DIR / "{lineage}_repeats.bed"    
    log:
        "logs/references/repeats/repeatmasker_combine_{lineage}.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        cat {input.simple} {input.complx} {input.known} {input.unknown} \
        | bedtools sort \
        | bedtools merge -c 4 -o collapse \
        | awk '{{print $1"\t"$2"\t"$3"\t"$4}}' > {output} 2> {log}
        """

# =================================================================================================
#   Per sample | Intercept depth by windows with repeats and call CNVs
# =================================================================================================

def cnv_calling_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "depth": INT_SAMPLES_DIR / "depth_quality" / s["sample"] / "depth_by_windows.tsv",
        "repeats": REFS_DIR / (s["lineage"] + "_repeats.bed")
        }
rule cnv_calling:
    input:
        unpack(cnv_calling_input)
    output:
        SAMPLES_DIR / "cnv" / "{sample}" / "cnv_calls.tsv"
    params:
        window_size = config["depth_quality"]["mosdepth"]["window"],
        depth_threshold = config["cnv"]["depth_threshold"]
    log:
        "logs/samples/cnv/cnv_calling_{sample}.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/cnv_calling.xsh "
        "-di {input.depth} "
        "-ri {input.repeats} "
        "-co {output} "
        "-sp {wildcards.sample} "
        "-wp {params.window_size} "
        "-dp {params.depth_threshold} &> {log}"
