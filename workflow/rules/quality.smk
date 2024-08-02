# =================================================================================================
#   Per sample | Run Mosdepth to get depth per region of all (raw) and good quality reads
# =================================================================================================

rule mosdepth:
    input:
        bam = OUTDIR / "snippy" / "{sample}" / "snps.bam",
        bai = OUTDIR / "snippy" / "{sample}" / "snps.bam.bai"
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
        bam = OUTDIR / "snippy" / "{sample}" / "snps.bam",
        bai = OUTDIR / "snippy" / "{sample}" / "snps.bam.bai"
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
        global_mode = OUTDIR / "depth_quality" / "{sample}" / "global_mode.tsv"
    output:
        OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/samples/depth_quality/depth_by_chrom_good_{sample}.log"
    shell:
        "python workflow/scripts/depth_by_chrom.py -b {input.bed} -g {input.global_mode} -o {output} -s {wildcards.sample} &> {log}"

# =================================================================================================
#   Per sample | Get mean MAPQ by region
# =================================================================================================

rule mapq:
    input:
       OUTDIR / "snippy" / "{sample}" / "snps.bam",
       rules.mosdepth.output.bed
    output:
        bed = OUTDIR / "depth_quality" / "{sample}" / "mapq.bed",
        region_bed = OUTDIR / "depth_quality" / "{sample}" / "mapq_region.bed"
    conda:
        "../envs/samtools.yaml"
    log:
       "logs/samples/depth_quality/mapq_{sample}.log"
    script:
        "../scripts/pileup_mapq.sh"