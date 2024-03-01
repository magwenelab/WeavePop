
def intersect_variants_input(wildcards):
    g = SAMPLE_REFERENCE[(SAMPLE_REFERENCE['group'] == wildcards.lineage)]
    gSAMPLES = list(g["sample"])
    return {
        "samples_vcfs": [OUTDIR / "snippy" / sample / "snps.vcf.gz" for sample in gSAMPLES],
    }

# Intersect the variants of all samples of the same lineage/group 
rule intersect_variants:
    input: 
        unpack(intersect_variants_input),
        gff = rules.gff2tsv.output,
    output: 
        DATASET_OUTDIR / "snps" / "{lineage}_variants.tsv",
    conda:
        "../envs/samtools.yaml"
    params:
        temp = DATASET_OUTDIR / "snps/temp_",
    log:
        "logs/snps/group_variants_{lineage}.log"
    shell:
        "xonsh workflow/scripts/variant_intersection.xsh -g {input.gff} -l {wildcards.lineage} -o {output} -t {params.temp}{wildcards.lineage} {input.samples_vcfs} &> {log}"