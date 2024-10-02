# =================================================================================================
#   Per lineage | Run RepeatModeler and RepeatMasker
# =================================================================================================


rule repeat_modeler_build:
    input:
        rules.links.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "RM_db" / "{lineage}.nsq",
    params:
        repdir=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats",
        name=lambda wildcards: INT_REFS_DIR
        / wildcards.lineage
        / "repeats"
        / "RM_db"
        / wildcards.lineage,
    log:
        "logs/references/repeats/repeatmodeler_build_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "mkdir -p {params.repdir}/RM_db && "
        "BuildDatabase "
        "-name {params.name} "
        "-engine ncbi "
        "{input} "
        "&> {log}"


rule repeat_modeler:
    input:
        database=rules.repeat_modeler_build.output,
        fasta=rules.links.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "RM_db" / "{lineage}-families.fa",
    params:
        repdir=lambda wildcards: INT_REFS_DIR
        / wildcards.lineage
        / "repeats"
        / "RModeler",
        db=lambda wildcards: INT_REFS_DIR
        / wildcards.lineage
        / "repeats"
        / "RM_db"
        / wildcards.lineage,
    log:
        "logs/references/repeats/repeatmodeler_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "RepeatModeler "
        "-database {params.db} "
        "-engine ncbi "
        "-pa {threads} "
        "-dir {params.repdir} "
        "&> {log}"


rule repeat_modeler_separate:
    input:
        fasta=rules.repeat_modeler.output,
    output:
        known=INT_REFS_DIR / "{lineage}" / "repeats" / "known.fasta",
        unknown=INT_REFS_DIR / "{lineage}" / "repeats" / "unknown.fasta",
    log:
        "logs/references/repeats/repeatmodeler_separate_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        cat {input} | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > {output.known} 2> {log}
        cat {input} | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > {output.unknown} 2>> {log}
        """


rule repeat_masker_1:
    input:
        database=config["cnv"]["repeats"]["repeats_database"],
        fasta=rules.links.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "01_simple" / "{lineage}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats" / "01_simple",
    log:
        "logs/references/repeats/repeatmasker1_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "RepeatMasker "
        "-pa {threads} "
        "-lib {input.database} "
        "-a "
        "-e ncbi "
        "-dir {params.dir} "
        "-noint "
        "-xsmall "
        "{input.fasta} "
        "&> {log}"


rule repeat_masker_2:
    input:
        database=config["cnv"]["repeats"]["repeats_database"],
        fasta=rules.links.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "02_complex" / "{lineage}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR
        / wildcards.lineage
        / "repeats"
        / "02_complex",
    log:
        "logs/references/repeats/repeatmasker2_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "RepeatMasker "
        "-pa {threads} "
        "-lib {input.database} "
        "-a "
        "-e ncbi "
        "-dir {params.dir} "
        "-nolow "
        "{input.fasta} "
        "&> {log}"


rule repeat_masker_3:
    input:
        known=rules.repeat_modeler_separate.output.known,
        fasta=rules.links.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "03_known" / "{lineage}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats" / "03_known",
    log:
        "logs/references/repeats/repeatmasker3_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "RepeatMasker "
        "-pa {threads} "
        "-lib {input.known} "
        "-a "
        "-e ncbi "
        "-dir {params.dir} "
        "-nolow {input.fasta} "
        "&> {log}"


rule repeat_masker_4:
    input:
        unknown=rules.repeat_modeler_separate.output.unknown,
        fasta=rules.links.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "04_unknown" / "{lineage}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR
        / wildcards.lineage
        / "repeats"
        / "04_unknown",
    log:
        "logs/references/repeats/repeatmasker4_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "RepeatMasker "
        "-pa {threads} "
        "-lib {input.unknown} "
        "-a "
        "-e ncbi "
        "-dir {params.dir} "
        "-nolow "
        "{input.fasta} "
        "&> {log}"


rule repeat_masker_bed:
    input:
        simple=rules.repeat_masker_1.output,
        complx=rules.repeat_masker_2.output,
        known=rules.repeat_masker_3.output,
        unknown=rules.repeat_masker_4.output,
    output:
        simple=INT_REFS_DIR / "{lineage}" / "repeats" / "01_simple" / "{lineage}.bed",
        complx=INT_REFS_DIR / "{lineage}" / "repeats" / "02_complex" / "{lineage}.bed",
        known=INT_REFS_DIR / "{lineage}" / "repeats" / "03_known" / "{lineage}.bed",
        unknown=INT_REFS_DIR / "{lineage}" / "repeats" / "04_unknown" / "{lineage}.bed",
    log:
        "logs/references/repeats/repeatmasker_combine_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        tail -n +4 {input.simple} | awk '{{print $5"\t"($6-1)"\t"$7"\t"$11}}' \
        1> {output.simple} 2> {log}
        tail -n +4 {input.complx} | awk '{{print $5"\t"($6-1)"\t"$7"\t"$11}}' \
        1> {output.complx} 2>> {log}
        tail -n +4 {input.known} | awk '{{print $5"\t"($6-1)"\t"$7"\t"$11}}' \
        1> {output.known} 2>> {log}
        tail -n +4 {input.unknown} | awk '{{print $5"\t"($6-1)"\t"$7"\t"$11}}' \
        1> {output.unknown} 2>> {log}
        """


rule repeat_masker_combine:
    input:
        simple=rules.repeat_masker_bed.output.simple,
        complx=rules.repeat_masker_bed.output.complx,
        known=rules.repeat_masker_bed.output.known,
        unknown=rules.repeat_masker_bed.output.unknown,
    output:
        REFS_DIR / "{lineage}_repeats.bed",
    log:
        "logs/references/repeats/repeatmasker_combine_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        cat {input.simple} {input.complx} {input.known} {input.unknown} \
        | bedtools sort \
        | bedtools merge -c 4 -o collapse \
        | awk '{{print $1"\t"$2"\t"$3"\t"$4}}' > {output} 2> {log}
        """
