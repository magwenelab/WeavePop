# =================================================================================================
#   Per sample | Run Mosdepth to get depth per region of all (raw) and good quality reads
# =================================================================================================

rule mosdepth:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai
    output:
        bed = OUTDIR / "mosdepth" / "{sample}" / "coverage.regions.bed.gz"
    params:
        window = config["depth_quality"]["mosdepth"]["window"],
        extra = config["depth_quality"]["mosdepth"]["extra"],
        outdir = OUTDIR / "mosdepth"
    conda: 
        "../envs/depth.yaml"
    threads:
      config["depth_quality"]["mosdepth"]["threads"]    
    log:
        "logs/samples/mosdepth/mosdepth_{sample}.log"
    shell:
        "mosdepth -n --by {params.window} -t {threads} {params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage {input.bam} "
        "&> {log}"

rule mosdepth_good:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai
    output:
        bed = OUTDIR / "mosdepth" / "{sample}" / "coverage_good.regions.bed.gz"
    params:
        window = config["depth_quality"]["mosdepth"]["window"],
        extra = config["depth_quality"]["mosdepth"]["extra"],
        min_mapq = config["depth_quality"]["mosdepth"]["min_mapq"],
        outdir = OUTDIR / "mosdepth"
    conda: 
        "../envs/depth.yaml"
    threads:
       config["depth_quality"]["mosdepth"]["threads"]   
    log:
        "logs/samples/mosdepth/mosdepth_good_{sample}.log"
    shell:
        "mosdepth -n --by {params.window} --mapq {params.min_mapq} -t {threads} {params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage_good {input.bam} "
        "&> {log}"

# =================================================================================================
#   Per sample | Get distribution and global mode fo depth (genome-wide depth to normalize)
# =================================================================================================

rule bam_good:
    input:
        bam = rules.snippy.output.bam
    output:
        bam_good = temp(OUTDIR / "depth_quality" / "{sample}" / "snps_good.bam"),
        bai_good = temp(OUTDIR / "depth_quality" / "{sample}" / "snps_good.bam.bai")
    conda:
        "../envs/samtools.yaml"
    params:
        min_mapq = config["depth_quality"]["mosdepth"]["min_mapq"]   
    log:
        "logs/samples/depth_quality/bam_good_{sample}.log"
    shell:
        "samtools view -q {params.min_mapq} -b {input} > {output.bam_good} 2> {log} && "
        "samtools index {output.bam_good} -o {output.bai_good} 2>> {log} "

def depth_distribution_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
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
        distrib = OUTDIR / "depth_quality" / "{sample}" / "depth_distribution.tsv",
        global_mode = temp(OUTDIR / "depth_quality" / "{sample}" / "global_mode.tsv")
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/samples/depth_quality/depth_distribution_{sample}.log"
    shell:
        "xonsh workflow/scripts/depth_distribution.xsh -s {wildcards.sample} -b {input.bam} -g {input.bam_good} -do {output.distrib} -go {output.global_mode} &> {log}"

# =================================================================================================
#   Per sample | Get depth by chromosome raw and good
# =================================================================================================

rule depth_by_chrom_raw:
    input:
        bed = rules.mosdepth.output.bed
    output:
        OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_raw.tsv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/samples/depth_quality/depth_by_chrom_raw_{sample}.log"
    shell:
        "python workflow/scripts/depth_by_chrom.py -b {input.bed} -o {output} -s {wildcards.sample} &> {log}"

rule depth_by_chrom_good:
    input:
        bed = rules.mosdepth_good.output.bed,
        global_mode = rules.depth_distribution.output.global_mode
    output:
        OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/samples/depth_quality/depth_by_chrom_good_{sample}.log"
    shell:
        "python workflow/scripts/depth_by_chrom.py -b {input.bed} -g {input.global_mode} -o {output} -s {wildcards.sample} &> {log}"



# =================================================================================================
#   Per sample | Get mapping stats
# =================================================================================================

rule mapping_stats:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai,
        global_mode = rules.depth_distribution.output.global_mode
    output:
        OUTDIR / "depth_quality" / "{sample}" / "mapping_stats.tsv"
    params:
        min_depth = config["depth_quality"]["flag_quality"]["min_percent_genome-wide_depth"],
        min_mapq = config["depth_quality"]["flag_quality"]["min_percent_MAPQ"],    
        min_pp= config["depth_quality"]["flag_quality"]["min_percent_properly_paired_reads"],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samples/depth_quality/mapping_stats_{sample}.log"
    shell:
        "xonsh workflow/scripts/mapping_stats.xsh -b {input.bam} -s {wildcards.sample} -m {input.global_mode} -d {params.min_depth} -q {params.min_mapq} -p {params.min_pp} -o {output} &> {log}"

# =================================================================================================
#   Per sample | Get mean MAPQ by region
# =================================================================================================

rule mapq:
    input:
       rules.snippy.output.bam,
       rules.mosdepth.output.bed
    output:
        bed = temp(OUTDIR / "depth_quality" / "{sample}" / "mapq.bed"),
        region_bed = temp(OUTDIR / "depth_quality" / "{sample}" / "mapq_region.bed")
    conda:
        "../envs/samtools.yaml"
    log:
       "logs/samples/depth_quality/mapq_{sample}.log"
    script:
        "../scripts/pileup_mapq.sh"