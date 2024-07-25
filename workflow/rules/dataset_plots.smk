# =================================================================================================
#   Per dataset | Join depth by chrom and mapping stats and plot summary
# =================================================================================================

rule join_depth_by_chrom_raw:
    input:
        expand(OUTDIR / "mosdepth" / "{sample}" / "depth_by_chrom_raw.tsv",sample=SAMPLES),
    output:
        DATASET_OUTDIR / "files" / "depth_by_chrom_raw.tsv",
    log:
        "logs/dataset/files/join_depth_by_chrom_raw.log"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp
        """

rule join_depth_by_chrom_good:
    input:
        expand(OUTDIR / "mosdepth" / "{sample}" / "depth_by_chrom_good.tsv",sample=SAMPLES),
    output:
        DATASET_OUTDIR / "files" / "depth_by_chrom_good.tsv",
    log:
        "logs/dataset/files/join_depth_by_chrom_good.log"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp
        """

rule join_mapping_stats:
    input:
        expand(OUTDIR / "samtools" / "{sample}" / "mapping_stats.tsv",sample=SAMPLES),
    output:
        DATASET_OUTDIR / "files" / "mapping_stats.tsv",
    log:
        "logs/dataset/files/join_mapping_stats.log"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp
        """

rule dataset_summary_plot:  
    input:
        SAMPLEFILE,
        CHROM_NAMES,
        rules.join_depth_by_chrom_good.output,
        rules.join_depth_by_chrom_raw.output,
        rules.join_mapping_stats.output
    output:
        DATASET_OUTDIR / "plots" / "dataset_summary.png"
    conda:
        "../envs/r.yaml"
    params:
        config["plotting"]["scale"]
    log:
        "logs/dataset/plots/dataset_summary.log"
    script:
        "../scripts/dataset_summary_plot.R"

# =================================================================================================
#   Per dataset | Join normalized depth by chrom and plot
# =================================================================================================
rule join_depth_by_chrom_normalized:
    input:
        expand(OUTDIR / "mosdepth" / "{sample}" / "depth_by_chrom_good_normalized.tsv",sample=SAMPLES),
    output:
        DATASET_OUTDIR / "files" / "depth_by_chrom_good_normalized.tsv",
    log:
        "logs/dataset/files/join_depth_by_chrom_good_normalized.log"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp
        """

rule dataset_depth_by_chrom_plot:
    input:
        SAMPLEFILE,
        CHROM_NAMES,
        rules.join_depth_by_chrom_normalized.output
    output:
        DATASET_OUTDIR / "plots" / "dataset_depth_by_chrom.png"
    conda:
        "../envs/r.yaml"
    params:
        column = config["plotting"]["metadata2color"],
        scale = config["plotting"]["scale"]
    log:
        "logs/dataset/plots/dataset_depth_by_chrom.log"
    script:
        "../scripts/dataset_depth_by_chrom_plot.R"