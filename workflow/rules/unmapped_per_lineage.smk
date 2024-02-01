rule unmapped_samples_plot:
    input:
        SAMPLEFILE,
        expand(rules.gff2tsv.output, lineage=LINEAGES),
        expand(rules.liftoff.output.unmapped, sample=SAMPLES)        
    output:
        DATASET_OUTDIR / "files" / "unmapped_count.tsv",
        DATASET_OUTDIR / "plots" / "unmapped.png"
    conda:
        "../envs/r.yaml"
    params:
        dir = REFDIR,
        sampledir = OUTDIR / "liftoff"
    log:
        "logs/liftoff/unmapped_count_plot.log"
    script:
        "../scripts/samples_unmapped_per_lin.R"