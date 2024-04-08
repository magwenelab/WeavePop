rule unmapped_samples_plot:
    input:
        SAMPLEFILE,
        expand(rules.agat_make_gff.output.tsv, lineage=LINEAGES),
        expand(rules.liftoff.output.unmapped, sample=SAMPLES)        
    output:
        expand(DATASET_OUTDIR / "files" / "{lineage}_unmapped_count.tsv", lineage=LINEAGES),
        expand(DATASET_OUTDIR / "plots" / "{lineage}_unmapped.svg", lineage=LINEAGES)
    conda:
        "../envs/r.yaml"
    params:
        dir = REFDIR,
        sampledir = OUTDIR / "liftoff",
        datasetplots = DATASET_OUTDIR / "plots",
        datasetfiles = DATASET_OUTDIR / "files"
    log:
        "logs/liftoff/unmapped_samples_plot.log"
    script:
        "../scripts/samples_unmapped_per_lin.R"