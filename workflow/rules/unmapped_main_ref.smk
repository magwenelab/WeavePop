# =================================================================================================
#   Join lineages | Plot and count the features that were not lifted over from the main reference
# =================================================================================================

rule unmapped_ref_features:
    input:
        INT_DATASET_DIR / "metadata.csv",
        rules.fix_gff_tsv.output.tsv,
        expand(rules.ref2ref_liftoff.output.unmapped, lineage=LINEAGES)        
    output:
        REFS_DIR / "lifotff" / "unmapped_count.tsv",
        REFS_DIR / "lifotff" / "unmapped.svg"
    conda:
        "../envs/r.yaml"
    params:
        refdir = REFS_DIR
    log:
        "logs/references/lifotff/unmapped_ref_features.log"
    script:
        "../scripts/unmapped_features_refs.R"

# =================================================================================================
#   Per dataset | Plot and count the features that were not lifted over
# =================================================================================================

rule unmapped_samples_plot:
    input:
        INT_DATASET_DIR / "metadata.csv",
        rules.fix_gff_tsv.output.tsv,
        expand(rules.liftoff.output.unmapped, sample=SAMPLES)        
    output:
        DATASET_DIR / "liftoff" / "unmapped_count.tsv",
        DATASET_DIR / "liftoff" / "unmapped.svg"
    conda:
        "../envs/r.yaml"
    params:
        dir = SAMPLES_DIR / "liftoff"
    log:
        "logs/dataset/liftoff/unmapped_samples_plot.log"
    script:
        "../scripts/samples_unmapped_main.R"