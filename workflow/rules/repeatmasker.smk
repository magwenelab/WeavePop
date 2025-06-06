# =================================================================================================
#   Per lineage | Run RepeatModeler and RepeatMasker
# =================================================================================================


rule repeat_modeler_build:
    input:
        rules.ref_fasta_symlinks.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "db_rmodeler" / "{lineage}.nsq",
    params:
        repdir=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats",
        name=lambda wildcards: INT_REFS_DIR
        / wildcards.lineage
        / "repeats"
        / "db_rmodeler"
        / wildcards.lineage,
    log:
        LOGS / "references" / "repeats" / "repeatmodeler_build_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "mkdir -p {params.repdir}/db_rmodeler && "
        "BuildDatabase "
        "-name {params.name} "
        "{input} "
        "&> {log}"


rule repeat_modeler:
    input:
        database=rules.repeat_modeler_build.output,
        fasta=rules.ref_fasta_symlinks.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "db_rmodeler" / "{lineage}-families.fa",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmodeler_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "wd=$(pwd) && "
        "cd $wd/{params.dir} && "
        "export HOME={params.dir} && "
        "RepeatModeler "
        "-database db_rmodeler/{wildcards.lineage} "
        "-engine ncbi "
        "-pa {threads} "
        "&> $wd/{log}"


rule repeat_modeler_separate:
    input:
        fasta=rules.repeat_modeler.output,
    output:
        known=INT_REFS_DIR / "{lineage}" / "repeats" / "known.fasta",
        unknown=INT_REFS_DIR / "{lineage}" / "repeats" / "unknown.fasta",
    log:
        LOGS / "references" / "repeats" / "repeatmodeler_separate_{lineage}.log",
    conda:
        "../envs/agat.yaml"
    shell:
        """
        cat {input} | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > {output.known} 2> {log} || true
        cat {input} | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > {output.unknown} 2>> {log} || true
        """


rule repeat_masker_1:
    input:
        database=REPEATS_FILE,
        fasta=rules.ref_fasta_symlinks.output,
        rule_order=rules.repeat_modeler_separate.output.known,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "01_simple" / "{lineage}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats" / "01_simple",
        tmp=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmasker1_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "wd=$(pwd) && "
        "cd {params.tmp} && "
        "export HOME={params.dir} && "
        "RepeatMasker "
        "-pa {threads} "
        "-lib $wd/{input.database} "
        "-a "
        "-e ncbi "
        "-dir $wd/{params.dir} "
        "-noint "
        "-xsmall "
        "$wd/{input.fasta} "
        "&> $wd/{log}"


rule repeat_masker_2:
    input:
        database=REPEATS_FILE,
        fasta=rules.ref_fasta_symlinks.output,
        rule_order=rules.repeat_masker_1.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "02_complex" / "{lineage}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats" / "02_complex",
        tmp=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmasker2_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "wd=$(pwd) && "
        "cd {params.tmp} && "
        "export HOME={params.dir} && "
        "RepeatMasker "
        "-pa {threads} "
        "-lib $wd/{input.database} "
        "-a "
        "-e ncbi "
        "-dir $wd/{params.dir} "
        "-nolow "
        "$wd/{input.fasta} "
        "&> $wd/{log}"


rule repeat_masker_3:
    input:
        known=rules.repeat_modeler_separate.output.known,
        fasta=rules.ref_fasta_symlinks.output,
        rule_order=rules.repeat_masker_2.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "03_known" / "{lineage}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats" / "03_known",
        tmp=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmasker3_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "mkdir -p $(dirname {output}) && "
        "if [ ! -s {input.known} ]; then "
        "touch {output} && "
        "echo 'Skipping RepeatMasker of known repeats because no known families were identified' "
        "&> {log}; "
        "else "
        "wd=$(pwd) && "
        "cd {params.tmp} && "
        "export HOME={params.dir} && "
        "RepeatMasker "
        "-pa {threads} "
        "-lib $wd/{input.known} "
        "-a "
        "-e ncbi "
        "-dir $wd/{params.dir} "
        "-nolow $wd/{input.fasta} "
        "&> $wd/{log}; "
        "fi"


rule repeat_masker_4:
    input:
        unknown=rules.repeat_modeler_separate.output.unknown,
        fasta=rules.ref_fasta_symlinks.output,
        rule_order=rules.repeat_masker_3.output,
    output:
        INT_REFS_DIR / "{lineage}" / "repeats" / "04_unknown" / "{lineage}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats" / "04_unknown",
        tmp=lambda wildcards: INT_REFS_DIR / wildcards.lineage / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmasker4_{lineage}.log",
    threads: config["cnv"]["repeats"]["repeats_threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "mkdir -p $(dirname {output}) && "
        "if [ ! -s {input.unknown} ]; then "
        "touch {output} && "
        "echo 'Skipping RepeatMasker of unknown repeats because no unknown families were identified' "
        "&> {log}; "
        "else "
        "wd=$(pwd) && "
        "cd {params.tmp} && "
        "export HOME={params.dir} && "
        "RepeatMasker "
        "-pa {threads} "
        "-lib $wd/{input.unknown} "
        "-a "
        "-e ncbi "
        "-dir $wd/{params.dir} "
        "-nolow "
        "$wd/{input.fasta} "
        "&> $wd/{log}; "
        "fi"


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
        LOGS / "references" / "repeats" / "repeatmasker_combine_{lineage}.log",
    conda:
        "../envs/shell.yaml"
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
        REFS_DIR / "{lineage}" / "{lineage}_repeats.bed",
    log:
        LOGS / "references" / "repeats" / "repeatmasker_combine_{lineage}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        cat {input.simple} {input.complx} {input.known} {input.unknown} \
        | bedtools sort \
        | bedtools merge -c 4 -o collapse \
        | awk '{{print $1"\t"$2"\t"$3"\t"$4}}' > {output} 2> {log}
        """
