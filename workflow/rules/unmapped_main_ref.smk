# Plot and count the features that were not lifted over from the main reference
rule unmapped_ref_features:
    input:
        SAMPLEFILE,
        rules.fix_gff_tsv.output.tsv,
        expand(rules.ref2ref_liftoff.output.unmapped, lineage=LINEAGES)        
    output:
        REFDIR / "unmapped_count.tsv",
        REFDIR / "unmapped.svg"
    conda:
        "../envs/r.yaml"
    params:
        refdir = REFDIR
    log:
        "logs/references/unmapped_ref_features.log"
    script:
        "../scripts/unmapped_features_refs.R"

rule unmapped_samples_plot:
    input:
        SAMPLEFILE,
        rules.fix_gff_tsv.output.tsv,
        expand(rules.liftoff.output.unmapped, sample=SAMPLES)        
    output:
        DATASET_OUTDIR / "files" / "unmapped_count.tsv",
        DATASET_OUTDIR / "plots" / "unmapped.svg"
    conda:
        "../envs/r.yaml"
    params:
        dir = OUTDIR / "liftoff"
    log:
        "logs/liftoff/unmapped_samples_plot.log"
    script:
        "../scripts/samples_unmapped_main.R"