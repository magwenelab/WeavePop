
def intersect_variants_input(wildcards):
    g = SAMPLE_REFERENCE[(SAMPLE_REFERENCE['group'] == wildcards.lineage)]
    gSAMPLES = list(g["sample"])
    return {
        "samples_vcfs": [OUTDIR / "snippy" / sample / "snps.vcf.gz" for sample in gSAMPLES]
    }

# Intersect the variants of all samples of the same lineage/group 
rule intersect_variants:
    input: 
        unpack(intersect_variants_input)
    output: 
        vcf = DATASET_OUTDIR / "snps" / "{lineage}_variants.vcf",
    conda:
        "../envs/samtools.yaml"
    params:
        sample_names=lambda wildcards, input: input.samples_vcfs
    log:
        "logs/snps/group_variants_{lineage}.log"
    shell:
        """
        echo {params.sample_names} > {output.sample_names} 2> {log} && 
        bcftools isec {input.samples_vcfs} 1> {output.vcf} 2>> {log}
        """
