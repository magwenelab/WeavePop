# =================================================================================================
#   Per sample | Get distribution and genome-wide depth to normalize
# =================================================================================================


rule bam_good:
    input:
        bam=rules.snippy.output.bam,
    output:
        bam_good=INT_SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "snps_good.bam",
        bai_good=INT_SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "snps_good.bam.bai",
    params:
        min_mapq=config["depth_quality"]["flag_quality"]["min_mapq"],
    log:
        LOGS / "samples" / "depth_quality" / "bam_good_{unf_sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -q {params.min_mapq} -b {input} > {output.bam_good} 2> {log} && "
        "samtools index {output.bam_good} -o {output.bai_good} 2>> {log} "


rule depth_distribution:
    input:
        unpack(depth_distribution_input),
    output:
        distrib=INT_SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "depth_distribution.tsv",
        summary=INT_SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "depth_summary.tsv",
    log:
        LOGS / "samples" / "depth_quality" / "depth_distribution_{unf_sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/depth_distribution.xsh "
        "-s {wildcards.unf_sample} "
        "-b {input.bam} "
        "-g {input.bam_good} "
        "-do {output.distrib} "
        "-so {output.summary} &> {log} "


# =================================================================================================
#   Per sample | Get mapping stats
# =================================================================================================


rule mapping_stats:
    input:
        bam=rules.snippy.output.bam,
        bai=rules.snippy.output.bai,
        depth=rules.depth_distribution.output.summary,
    output:
        SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "mapping_stats.tsv",
    params:
        low_mapq=config["depth_quality"]["flag_quality"]["low_MAPQ_limit"],
        high_mapq=config["depth_quality"]["flag_quality"]["high_MAPQ_limit"],
        min_mapq=config["depth_quality"]["flag_quality"]["min_mapq"],
    log:
        LOGS / "samples" / "depth_quality" / "mapping_stats_{unf_sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/mapping_stats.xsh "
        "-s {wildcards.unf_sample} "
        "-b {input.bam} "
        "-d {input.depth} "
        "-l {params.low_mapq} "
        "-h {params.high_mapq} "
        "-mq {params.min_mapq} "
        "-o {output} &> {log}"


# =================================================================================================
#   Per dataset | Join mapping stats
# =================================================================================================


rule join_mapping_stats:
    input:
        expand(
            SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "mapping_stats.tsv",
            unf_sample=UNFILT_SAMPLES,
        ),
    output:
        DATASET_DIR / "depth_quality" / "mapping_stats.tsv",
    params:
        min_depth=config["depth_quality"]["flag_quality"]["min_percent_genome-wide_depth"],
        min_high_mapq=config["depth_quality"]["flag_quality"]["min_percent_MAPQ"],
        min_pm=config["depth_quality"]["flag_quality"]["min_percent_mapped_reads"],
        min_coverage=config["depth_quality"]["flag_quality"]["min_percent_coverage"],
    log:
        LOGS / "samples" / "depth_quality" / "join_mapping_stats.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_mapping_stats.py"

# =================================================================================================
#   All refernces | Obtain chromosome lengths
# =================================================================================================


rule chromosome_lengths:
    input:
        INT_REFS_DIR / "{unf_lineage}" / "{unf_lineage}.fasta",         
    output:
        INT_REFS_DIR / "{unf_lineage}" / "chromosome_lengths.tsv",
    log:
        LOGS / "references" / "depth_quality" / "{unf_lineage}_chromosome_lengths.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        """
        seqkit fx2tab -l -i -n {input} |\
        awk -v lin={wildcards.unf_lineage} '{{print lin, $0}}' OFS='\t' \
        1> {output} 2> {log}
        """

rule join_chromosome_lengths:
    input:
        chrom_names=CHROM_NAMES,  
        chrom_lengths=expand(INT_REFS_DIR / "{unf_lineage}" / "chromosome_lengths.tsv", unf_lineage=UNF_LINEAGES),
    output:
        INT_REFS_DIR / "chromosome_lengths.tsv",
    log:
        LOGS / "references" / "ref_processing" / "join_chromosome_lengths.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_chromosome_lengths.py"        

# =================================================================================================
#   Per dataset | Checkpoint to filter out low quality samples
# =================================================================================================


rule quality_filter:
    input:
        rules.join_mapping_stats.output,
        rules.join_chromosome_lengths.output,
    output:
        metadata=DATASET_DIR / "metadata.csv",
        chromosomes=DATASET_DIR / "chromosomes.csv",
    params:
        filter=config["depth_quality"]["flag_quality"]["filter"],
        metadata=UNFILT_SAMPLE_TABLE
    log:
        LOGS / "samples" / "depth_quality" / "quality_filter.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/quality_filter.py"


checkpoint filter_wildcards:
    input:
        rules.quality_filter.output.metadata,
    output:
        directory(INT_SAMPLES_DIR / "filtered_samples"),
        directory(INT_REFS_DIR / "filtered_lineages"),
    log:
        LOGS / "samples" / "depth_quality" / "filter_wildcards.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/filter_wildcards.py"


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