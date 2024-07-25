# =================================================================================================
#   Join samples per lineage | Intersect VCF files of all samples from each lineage and annotate
# =================================================================================================
def intersect_vcfs_input(wildcards):
    l = LINEAGE_REFERENCE.loc[wildcards.lineage,]
    return {
        "vcfs" : expand(OUTDIR / "snippy" / "{sample}" / "snps.vcf.gz", sample=l["sample"])
    }
    
rule intersect_vcfs:
    input:
        unpack(intersect_vcfs_input)
    output:
        vcf = DATASET_OUTDIR / "snps" / "{lineage}_intersection.vcf",
        tsv = DATASET_OUTDIR / "snps" / "{lineage}_presence.tsv"
    params:
        tmp_dir = "{lineage}"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/intersect_vcfs_{lineage}.log"
    shell:
        "xonsh workflow/scripts/intersect_vcfs.xsh "
        "-v {output.vcf} "
        "-p {output.tsv} "
        "-l {wildcards.lineage} "
        "-t {params.tmp_dir} "
        "{input.vcfs} "
        "&> {log}"

rule snpeff:
    input:
        vcf = rules.intersect_vcfs.output.vcf,
        db_done = rules.build_refs_db.output
    output:
        vcf = DATASET_OUTDIR / "snps" / "{lineage}_snpeff.vcf",
        html = DATASET_OUTDIR / "snps" / "{lineage}_snpeff.html"
    params:
       dir = os.getcwd() / REFDIR / "snpeff_data",
       config = REFDIR / "snpeff_data" / "snpEff.config",
       name = config["species_name"] + "_{lineage}"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/snpeff_{lineage}.log"
    shell:
        "snpEff ann -v -classic "
        "-dataDir {params.dir} "
        "-config {params.config} "
        "-s {output.html} "
        "{params.name} "
        "{input.vcf} "
        "1> {output.vcf} 2> {log}"

rule extract_vcf_annotation:
    input:
        vcf = rules.snpeff.output.vcf
    output:
        effects = DATASET_OUTDIR / "snps" / "{lineage}_effects.tsv",
        variants = DATASET_OUTDIR / "snps" / "{lineage}_variants.tsv",
        lofs = DATASET_OUTDIR / "snps" / "{lineage}_lofs.tsv",
        nmds = DATASET_OUTDIR / "snps" / "{lineage}_nmds.tsv"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/extract_vcf_annotation_{lineage}.log"
    shell:
        "xonsh workflow/scripts/extract_vcf_annotation.xsh "
        "-i {input.vcf} "
        "-e {output.effects} "
        "-v {output.variants} "
        "-f {output.lofs} "
        "-n {output.nmds} "
        "-l {wildcards.lineage} "
        "&> {log}"

