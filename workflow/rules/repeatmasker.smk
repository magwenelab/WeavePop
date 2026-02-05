# =================================================================================================
#   Per ref_genome | Run RepeatModeler and RepeatMasker
# =================================================================================================


rule repeat_modeler_build:
    input:
        rules.ref_fasta_symlinks.output,
    output:
        INT_REFS_DIR / "{ref_genome}" / "repeats" / "db_rmodeler" / "{ref_genome}.nsq",
    params:
        repdir=lambda wildcards: INT_REFS_DIR / wildcards.ref_genome / "repeats",
        name=lambda wildcards: INT_REFS_DIR
        / wildcards.ref_genome
        / "repeats"
        / "db_rmodeler"
        / wildcards.ref_genome,
    log:
        LOGS / "references" / "repeats" / "repeatmodeler_build_{ref_genome}.log",
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
        INT_REFS_DIR / "{ref_genome}" / "repeats" / "db_rmodeler" / "{ref_genome}-families.fa",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.ref_genome / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmodeler_{ref_genome}.log",
    threads: config["repeats"]["threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        "wd=$(pwd) && "
        "cd $wd/{params.dir} && "
        "export HOME={params.dir} && "
        "RepeatModeler "
        "-database db_rmodeler/{wildcards.ref_genome} "
        "-engine ncbi "
        "-pa {threads} "
        "&> $wd/{log}"


rule repeat_modeler_separate:
    input:
        fasta=rules.repeat_modeler.output,
    output:
        known=INT_REFS_DIR / "{ref_genome}" / "repeats" / "known.fasta",
        unknown=INT_REFS_DIR / "{ref_genome}" / "repeats" / "unknown.fasta",
    log:
        LOGS / "references" / "repeats" / "repeatmodeler_separate_{ref_genome}.log",
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
        INT_REFS_DIR / "{ref_genome}" / "repeats" / "01_simple" / "{ref_genome}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.ref_genome / "repeats" / "01_simple",
        tmp=lambda wildcards: INT_REFS_DIR / wildcards.ref_genome / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmasker1_{ref_genome}.log",
    threads: config["repeats"]["threads"]
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
        INT_REFS_DIR / "{ref_genome}" / "repeats" / "02_complex" / "{ref_genome}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR
        / wildcards.ref_genome
        / "repeats"
        / "02_complex",
        tmp=lambda wildcards: INT_REFS_DIR / wildcards.ref_genome / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmasker2_{ref_genome}.log",
    threads: config["repeats"]["threads"]
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
        INT_REFS_DIR / "{ref_genome}" / "repeats" / "03_known" / "{ref_genome}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR / wildcards.ref_genome / "repeats" / "03_known",
        tmp=lambda wildcards: INT_REFS_DIR / wildcards.ref_genome / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmasker3_{ref_genome}.log",
    threads: config["repeats"]["threads"]
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
        INT_REFS_DIR / "{ref_genome}" / "repeats" / "04_unknown" / "{ref_genome}.fasta.out",
    params:
        dir=lambda wildcards: INT_REFS_DIR
        / wildcards.ref_genome
        / "repeats"
        / "04_unknown",
        tmp=lambda wildcards: INT_REFS_DIR / wildcards.ref_genome / "repeats",
    log:
        LOGS / "references" / "repeats" / "repeatmasker4_{ref_genome}.log",
    threads: config["repeats"]["threads"]
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
        simple=INT_REFS_DIR / "{ref_genome}" / "repeats" / "01_simple" / "{ref_genome}.bed",
        complx=INT_REFS_DIR / "{ref_genome}" / "repeats" / "02_complex" / "{ref_genome}.bed",
        known=INT_REFS_DIR / "{ref_genome}" / "repeats" / "03_known" / "{ref_genome}.bed",
        unknown=INT_REFS_DIR / "{ref_genome}" / "repeats" / "04_unknown" / "{ref_genome}.bed",
    log:
        LOGS / "references" / "repeats" / "repeatmasker_combine_{ref_genome}.log",
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
        bed=REFS_DIR / "{ref_genome}" / "{ref_genome}_repeats.bed",
    log:
        LOGS / "references" / "repeats" / "repeatmasker_combine_{ref_genome}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        cat {input.simple} {input.complx} {input.known} {input.unknown} \
        | bedtools sort \
        | bedtools merge -c 4 -o distinct \
        | awk '{{print $1"\t"$2"\t"$3"\t"$4}}' > {output.bed} 2> {log}
        """
