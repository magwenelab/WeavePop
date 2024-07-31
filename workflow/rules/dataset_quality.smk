# =================================================================================================
#   Per dataset | Join depth by chrom and mapping stats 
# =================================================================================================

rule join_depth_by_chrom_raw:
    input:
        expand(OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_raw.tsv",sample=SAMPLES),
    output:
        DATASET_OUTDIR / "depth_quality" / "depth_by_chrom_raw.tsv",
    log:
        "logs/dataset/depth_quality/join_depth_by_chrom_raw.log"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp
        """

rule join_depth_by_chrom_good:
    input:
        expand(OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv",sample=SAMPLES),
    output:
        DATASET_OUTDIR / "depth_quality" / "depth_by_chrom_good.tsv",
    log:
        "logs/dataset/depth_quality/join_depth_by_chrom_good.log"
    shell:
        """
        head -q -n 1 {input} 1> {output}.temp 2>> {log}
        head -n 1 {output}.temp 1> {output} 2>> {log}
        tail -q -n +2 {input} 1>> {output} 2>> {log}
        rm {output}.temp
        """

rule join_mapping_stats:
    input:
        expand(OUTDIR / "depth_quality" / "{sample}" / "mapping_stats.tsv",sample=SAMPLES),
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
