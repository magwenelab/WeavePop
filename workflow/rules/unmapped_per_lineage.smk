rule unmapped_samples_plot:
    input:
        SAMPLEFILE,
        expand(rules.gff2tsv.output),
        expand(rules.liftoff.output.unmapped, sample=SAMPLES)        
    output:
        DATASET_OUTDIR / "files" / "unmapped_count.tsv",
        DATASET_OUTDIR / "plots" / "unmapped.svg"
    conda:
        "../envs/r.yaml"
    params:
        script = workflow.source_path("../scripts/samples_unmapped_per_lin.R"),
        dir = REFDIR,
        sampledir = OUTDIR / "liftoff"
    log:
        "logs/liftoff/unmapped_count_plot.log"
    script:
        "{params.script}"