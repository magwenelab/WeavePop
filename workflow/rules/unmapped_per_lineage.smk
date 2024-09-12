# =================================================================================================
#   Per dataset | Plot and count the features that were not lifted over
# =================================================================================================

rule unmapped_samples_plot:
    input:
        INT_DATASET_DIR / "metadata.csv",
        expand(rules.fix_gff_tsv.output.tsv, lineage=LINEAGES),
        expand(rules.liftoff.output.unmapped, sample=SAMPLES)        
    output:
        expand(DATASET_DIR / "liftoff" / "{lineage}_unmapped_count.tsv", lineage=LINEAGES),
        expand(DATASET_DIR / "plots" / "{lineage}_unmapped.svg", lineage=LINEAGES)
    conda:
        "../envs/r.yaml"
    params:
        dir = REFS_DIR,
        sampledir = SAMPLES_DIR / "liftoff",
        datasetplots = DATASET_DIR / "plots",
        datasetfiles = DATASET_DIR / "liftoff"
    log:
        "logs/datset/plots/unmapped_samples_plot.log"
    script:
        "../scripts/samples_unmapped_per_lin.R"