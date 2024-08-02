# =================================================================================================
#   Join lineages | Plot and count the features that were not lifted over from the main reference
# =================================================================================================

rule unmapped_ref_features:
    input:
        GENERAL_OUTPUT / "metadata.csv",
        rules.fix_gff_tsv.output.tsv,
        expand(rules.ref2ref_liftoff.output.unmapped, lineage=LINEAGES)        
    output:
        REFDIR / "lifotff" / "unmapped_count.tsv",
        REFDIR / "lifotff" / "unmapped.svg"
    conda:
        "../envs/r.yaml"
    params:
        refdir = REFDIR
    log:
        "logs/references/lifotff/unmapped_ref_features.log"
    script:
        "../scripts/unmapped_features_refs.R"

# =================================================================================================
#   Per dataset | Plot and count the features that were not lifted over
# =================================================================================================

rule unmapped_samples_plot:
    input:
        GENERAL_OUTPUT / "metadata.csv",
        rules.fix_gff_tsv.output.tsv,
        expand(rules.liftoff.output.unmapped, sample=SAMPLES)        
    output:
        DATASET_OUTDIR / "liftoff" / "unmapped_count.tsv",
        DATASET_OUTDIR / "plots" / "unmapped.svg"
    conda:
        "../envs/r.yaml"
    params:
        dir = OUTDIR / "liftoff"
    log:
        "logs/dataset/plots/unmapped_samples_plot.log"
    script:
        "../scripts/samples_unmapped_main.R"