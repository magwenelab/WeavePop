# =================================================================================================
#   Per sample | Get distribution and global mode fo depth (genome-wide depth to normalize)
# =================================================================================================

rule bam_good:
    input:
        bam = rules.snippy.output.bam
    output:
        bam_good = INT_SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "snps_good.bam",
        bai_good = INT_SAMPLES_DIR/ "depth_quality" / "{unf_sample}" / "snps_good.bam.bai"
    params:
        min_mapq = config["depth_quality"]["mosdepth"]["min_mapq"]   
    log:
        "logs/samples/depth_quality/bam_good_{unf_sample}.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -q {params.min_mapq} -b {input} > {output.bam_good} 2> {log} && "
        "samtools index {output.bam_good} -o {output.bai_good} 2>> {log} "

rule depth_distribution:
    input:
        unpack(depth_distribution_input)
    output:
        distrib = INT_SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "depth_distribution.tsv",
        by_chrom_good = SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "depth_by_chrom_good.tsv",
        by_chrom_raw = SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "depth_by_chrom_raw.tsv"
    log:
        "logs/samples/depth_quality/depth_distribution_{unf_sample}.log"
    resources:
        tmpdir = TEMPDIR
    conda: 
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/depth_distribution.xsh "
        "-s {wildcards.unf_sample} "
        "-b {input.bam} "
        "-g {input.bam_good} " 
        "-do {output.distrib} " 
        "-go {output.by_chrom_good} " 
        "-ro {output.by_chrom_raw} &> {log} "

# =================================================================================================
#   Per sample | Get mapping stats
# =================================================================================================

rule mapping_stats:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai,
        global_mode = rules.depth_distribution.output.by_chrom_good
    output:
        SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "mapping_stats.tsv"
    params:
        low_mapq = config["depth_quality"]["flag_quality"]["low_MAPQ_limit"],
        high_mapq = config["depth_quality"]["flag_quality"]["high_MAPQ_limit"],
        min_position_depth = config["depth_quality"]["flag_quality"]["min_position_depth"],
        min_depth = config["depth_quality"]["flag_quality"]["min_percent_genome-wide_depth"],
        min_mapq = config["depth_quality"]["flag_quality"]["min_percent_MAPQ"],    
        min_pp= config["depth_quality"]["flag_quality"]["min_percent_properly_paired_reads"],
        min_coverage = config["depth_quality"]["flag_quality"]["min_percent_coverage"]
    log:
        "logs/samples/depth_quality/mapping_stats_{unf_sample}.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/mapping_stats.xsh "
        "-b {input.bam} "
        "-s {wildcards.unf_sample} "
        "-m {input.global_mode} "
        "-l {params.low_mapq} "
        "-h {params.high_mapq} "
        "-pd {params.min_position_depth} "
        "-d {params.min_depth} "
        "-q {params.min_mapq} "
        "-p {params.min_pp} "
        "-c {params.min_coverage} "
        "-o {output} &> {log}"

# =================================================================================================
#   Per dataset | Join mapping stats 
# =================================================================================================

rule join_mapping_stats:
    input:
        expand(SAMPLES_DIR / "depth_quality" / "{unf_sample}" / "mapping_stats.tsv",unf_sample=UNFILTERED_SAMPLES),
    output:
        INT_DATASET_DIR / "depth_quality" / "unfiltered_mapping_stats.tsv",
    log:
        "logs/dataset/depth_quality/join_mapping_stats.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/shell.yaml"
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
        stats = DATASET_DIR / "depth_quality" / "mapping_stats.tsv",
        metadata = INT_DATASET_DIR / "metadata.csv"
    params:
        exclude = config["depth_quality"]["flag_quality"]["exclude_samples"]
    log:
        "logs/dataset/depth_quality/quality_filter.log"
    run:
        with open(log[0], "w") as log_file:
            try:
                stats = pd.read_csv(input[0], sep="\t", header = 0)
                metadata = pd.read_csv(input[1], header=0)
                if params.exclude:
                    stats_filtered = stats[stats["quality_warning"].isna()]
                    metadata_filtered = metadata.loc[metadata["sample"].isin(stats_filtered["sample"]),]
                else:
                    stats_filtered = stats
                    metadata_filtered = metadata
                stats_filtered.to_csv(output.stats, index=False, header=True, sep = "\t")
                metadata_filtered.to_csv(output.metadata, index=False)
                log_file.write("Successfully filtered samples from tables.\n")
            except Exception as e:
                log_file.write(f"Error: {e}\n")
                raise e

checkpoint filtered_samples:
    input:
        rules.quality_filter.output.metadata
    output:
        directory(INT_SAMPLES_DIR / "filtered_samples")
    log:
        "logs/dataset/depth_quality/filtered_samples.log"
    run: 
        with open(log[0], "w") as log_file:
            try:
                metadata = pd.read_csv(input[0], header=0)
                sample_names = list(metadata["sample"])
                for sample_name in sample_names:
                    path_s = INT_SAMPLES_DIR / "filtered_samples" / f"{sample_name}.txt"
                    path_s.parent.mkdir(parents=True, exist_ok=True)
                    path_s.touch()
                log_file.write("Successfully created sample files.\n")
            except Exception as e:
                log_file.write(f"Error: {e}\n")
                raise e

checkpoint filtered_lineages:
    input:
        rules.quality_filter.output.metadata
    output:
        directory(INT_REFS_DIR / "filtered_lineages")
    log:
        "logs/dataset/depth_quality/filtered_lineages.log"
    run: 
        with open(log[0], "w") as log_file:
            try:
                lineages = list(pd.read_csv(input[0], header=0)["lineage"])
                for lineage in lineages:
                    path_l = INT_REFS_DIR / "filtered_lineages" / f"{lineage}.txt"
                    path_l.parent.mkdir(parents=True, exist_ok=True)
                    path_l.touch()
                log_file.write("Successfully created lineage files.\n")
            except Exception as e:
                log_file.write(f"Error: {e}\n")
                raise e

# =================================================================================================
#   Per dataset | Join depth by chrom 
# =================================================================================================

rule join_depth_by_chrom_raw:
    input:
        expand(SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_chrom_raw.tsv",sample=SAMPLES),
    output:
        DATASET_DIR / "depth_quality" / "depth_by_chrom_raw.tsv",
    log:
        "logs/dataset/depth_quality/join_depth_by_chrom_raw.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp 2>> {log}
        """

rule join_depth_by_chrom_good:
    input:
        expand(SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv",sample=SAMPLES),
    output:
        DATASET_DIR / "depth_quality" / "depth_by_chrom_good.tsv",
    log:
        "logs/dataset/depth_quality/join_depth_by_chrom_good.log"
    conda:
        "../envs/shell.yaml"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp 2>> {log}
        """
