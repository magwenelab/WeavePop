# =================================================================================================
#   Per sample | Get distribution and global mode fo depth (genome-wide depth to normalize)
# =================================================================================================

rule bam_good:
    input:
        bam = rules.snippy.output.bam
    output:
        bam_good = OUTDIR / "depth_quality" / "{unf_sample}" / "snps_good.bam",
        bai_good = OUTDIR / "depth_quality" / "{unf_sample}" / "snps_good.bam.bai"
    conda:
        "../envs/samtools.yaml"
    params:
        min_mapq = config["depth_quality"]["mosdepth"]["min_mapq"]   
    log:
        "logs/samples/depth_quality/bam_good_{unf_sample}.log"
    shell:
        "samtools view -q {params.min_mapq} -b {input} > {output.bam_good} 2> {log} && "
        "samtools index {output.bam_good} -o {output.bai_good} 2>> {log} "

def depth_distribution_input(wildcards):
    s = UNFILTERED_SAMPLE_REFERENCE.loc[wildcards.unf_sample,]
    return {
        "bam": OUTDIR / "snippy" / s["sample"] / "snps.bam" ,
        "bai": OUTDIR / "snippy" / s["sample"] / "snps.bam.bai",
        "bam_good": OUTDIR / "depth_quality" / s["sample"] / "snps_good.bam",
        "bai_good": OUTDIR / "depth_quality" / s["sample"] / "snps_good.bam.bai"
        }
rule depth_distribution:
    input:
        unpack(depth_distribution_input)
    output:
        distrib = OUTDIR / "depth_quality" / "{unf_sample}" / "depth_distribution.tsv",
        global_mode = OUTDIR / "depth_quality" / "{unf_sample}" / "global_mode.tsv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/samples/depth_quality/depth_distribution_{unf_sample}.log"
    shell:
        "xonsh workflow/scripts/depth_distribution.xsh -s {wildcards.unf_sample} -b {input.bam} -g {input.bam_good} -do {output.distrib} -go {output.global_mode} &> {log}"

# =================================================================================================
#   Per sample | Get mapping stats
# =================================================================================================

rule mapping_stats:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai,
        global_mode = rules.depth_distribution.output.global_mode
    output:
        OUTDIR / "depth_quality" / "{unf_sample}" / "mapping_stats.tsv"
    params:
        low_mapq = config["depth_quality"]["flag_quality"]["low_MAPQ_threshold"],
        high_mapq = config["depth_quality"]["flag_quality"]["high_MAPQ_threshold"],
        min_depth = config["depth_quality"]["flag_quality"]["min_percent_genome-wide_depth"],
        min_mapq = config["depth_quality"]["flag_quality"]["min_percent_MAPQ"],    
        min_pp= config["depth_quality"]["flag_quality"]["min_percent_properly_paired_reads"],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samples/depth_quality/mapping_stats_{unf_sample}.log"
    shell:
        "xonsh workflow/scripts/mapping_stats.xsh -b {input.bam} -s {wildcards.unf_sample} -m {input.global_mode} -l {params.low_mapq} -h {params.high_mapq} -d {params.min_depth} -q {params.min_mapq} -p {params.min_pp} -o {output} &> {log}"

# =================================================================================================
#   Per dataset | Join mapping stats 
# =================================================================================================

rule join_mapping_stats:
    input:
        expand(OUTDIR / "depth_quality" / "{unf_sample}" / "mapping_stats.tsv",unf_sample=UNFILTERED_SAMPLES),
    output:
        DATASET_OUTDIR / "depth_quality" / "mapping_stats.tsv",
    log:
        "logs/dataset/depth_quality/join_mapping_stats.log"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp
        """

# =================================================================================================
#   Per dataset | Checkpoint to filter out low quality samples
# =================================================================================================

rule quality_filter:
    input:
        rules.join_mapping_stats.output,
        UNFILTERED_SAMPLE_FILE
    output:
        stats = DATASET_OUTDIR / "depth_quality" / "filtered_mapping_stats.tsv",
        metadata = GENERAL_OUTPUT / "metadata.csv"
    run:
        stats = pd.read_csv(input[0], sep="\t", header = 0)
        stats_filtered = stats[stats["quality_warning"].isna()]
        stats_filtered.to_csv(output.stats, index=False, header=True, sep = "\t")
        metadata_unfiltered = pd.read_csv(input[1], header=0)
        metadata_filtered = metadata_unfiltered.loc[metadata_unfiltered["sample"].isin(stats_filtered["sample"]),]
        metadata_filtered.to_csv(output.metadata, index=False)

checkpoint filtered_samples:
    input:
        rules.quality_filter.output.metadata
    output:
        directory(GENERAL_OUTPUT / "filtered_samples")
    run: 
        metadata = pd.read_csv(input[0], header=0)
        sample_names = list(metadata["sample"])
        for sample_name in sample_names:
            path_s = GENERAL_OUTPUT / "filtered_samples" / f"{sample_name}.txt"
            path_s.parent.mkdir(parents=True, exist_ok=True)
            path_s.touch()
        
checkpoint filtered_lineages:
    input:
        rules.quality_filter.output.metadata
    output:
        directory(GENERAL_OUTPUT / "filtered_lineages")
    run: 
        lineages = list(pd.read_csv(input[0], header=0)["lineage"])
        for lineage in lineages:
            path_l = GENERAL_OUTPUT / "filtered_lineages" / f"{lineage}.txt"
            path_l.parent.mkdir(parents=True, exist_ok=True)
            path_l.touch()