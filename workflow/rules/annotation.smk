# =================================================================================================
# Per sample | Run Liftoff to annotate the assembly with the corresponding reference genome
# =================================================================================================


rule liftoff:
    input:
        unpack(liftoff_input),
    output:
        ref_gff=INT_SAMPLES_DIR / "annotation" / "{sample}" / "liftoff" / "ref.gff",
        ref_gff_db=temp(
            INT_SAMPLES_DIR / "annotation" / "{sample}" / "liftoff" / "ref.gff_db"
        ),
        gff=INT_SAMPLES_DIR / "annotation" / "{sample}" / "liftoff" / "lifted.gff",
        polished=INT_SAMPLES_DIR
        / "annotation"
        / "{sample}"
        / "liftoff"
        / "lifted.gff_polished",
        unmapped=INT_SAMPLES_DIR
        / "annotation"
        / "{sample}"
        / "liftoff"
        / "unmapped_features.txt",
        intermediate=temp(
            directory(
                INT_SAMPLES_DIR
                / "annotation"
                / "{sample}"
                / "liftoff"
                / "intermediate_liftoff"
            )
        ),
        fai=temp(
            SAMPLES_DIR / "mapping_and_variants" / "{sample}" / "snps.consensus.fa.fai"
        ),
        mmi=temp(
            SAMPLES_DIR / "mapping_and_variants" / "{sample}" / "snps.consensus.fa.mmi"
        ),
    params:
        extra=config["annotation"]["liftoff"]["extra"],
        outpath=INT_SAMPLES_DIR / "annotation" / "{sample}" / "liftoff",
    log:
        LOGS / "samples" / "annotation" / "liftoff_{sample}.log",
    threads: config["annotation"]["liftoff"]["threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/liftoff.yaml"
    shell:
        "ln -s -r -f {input.refgff} {output.ref_gff} &> {log} || true "
        "&& "
        "liftoff "
        "-g {output.ref_gff} "
        "-o {output.gff} "
        "-dir {output.intermediate} "
        "-u {output.unmapped} "
        "-p {threads} "
        "-polish "
        "{params.extra} "
        "{input.target} "
        "{input.refgenome} &>> {log}"


# =================================================================================================
# Per sample | Annotate intergenic regions and introns
# =================================================================================================


rule add_intergenic:
    input:
        gff=rules.liftoff.output.polished,
        config=rules.agat_config.output,
    output:
        gff=temp(INT_SAMPLES_DIR / "annotation" / "{sample}" / "intergenic.gff"),
    log:
        LOGS / "samples" / "annotation" / "add_intergenic_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_add_intergenic_regions.pl "
        "-g {input.gff} "
        "-o {output.gff} "
        "-c {input.config} "
        "&> {log} "


rule add_introns:
    input:
        gff=rules.add_intergenic.output.gff,
        config=rules.agat_config.output,
    output:
        gff=temp(INT_SAMPLES_DIR / "annotation" / "{sample}" / "interg_introns.gff"),
    log:
        LOGS / "samples" / "annotation" / "add_introns_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_add_introns.pl "
        "-g {input.gff} "
        "-o {output.gff} "
        "-c {input.config} "
        "&> {log} "


rule annotation_gff2tsv:
    input:
        gff=rules.add_introns.output.gff,
        config=rules.agat_config.output,
    output:
        tsv=temp(INT_SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff.tsv"),
    log:
        LOGS / "samples" / "annotation" / "annotation_gff2tsv_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "-g {input.gff} "
        "-o {output.tsv} "
        "-c {input.config} "
        "&> {log} "


rule reformat_annotation:
    input:
        tsv=rules.annotation_gff2tsv.output.tsv,
    output:
        tsv=SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff.tsv",
        gff=SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff",
    params:
        version="sample",
    log:
        LOGS / "samples" / "annotation" / "reformat_annotation_{sample}.log",
    conda:
        "../envs/shell.yaml"
    script:
        "../scripts/reformat_annotation.py"


# =================================================================================================
# Per sample | Run AGAT to extract CDS and protein sequences
# =================================================================================================


rule extract_cds:
    input:
        gff=rules.reformat_annotation.output.gff,
        fa=SAMPLES_DIR / "mapping_and_variants" / "{sample}" / "snps.consensus.fa",
        config=rules.agat_config.output,
    output:
        fa=SAMPLES_DIR / "annotation" / "{sample}" / "cds.fa",
    log:
        LOGS / "samples" / "annotation" / "extract_cds_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} "
        "-f {input.fa} "
        "-o {output.fa} "
        "-c {input.config} "
        "&> {log} "


rule extract_prots:
    input:
        gff=rules.reformat_annotation.output.gff,
        fa=SAMPLES_DIR / "mapping_and_variants" / "{sample}" / "snps.consensus.fa",
        config=rules.agat_config.output,
        cds=rules.extract_cds.output.fa,
    output:
        fa=SAMPLES_DIR / "annotation" / "{sample}" / "proteins.fa",
    log:
        LOGS / "samples" / "annotation" / "extract_prots_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} "
        "-f {input.fa} "
        "-o {output.fa} "
        "-p "
        "-c {input.config} &> {log} "
        "&& "
        "sed -i 's/type=cds//g' {output} &>> {log} "


# =================================================================================================
# Per sample | Convert fasta to csv to include in database
# =================================================================================================


rule cds2csv:
    input:
        fa=rules.extract_cds.output.fa,
    output:
        csv=INT_SAMPLES_DIR / "annotation" / "{sample}" / "cds.csv",
    log:
        LOGS / "samples" / "annotation" / "cds2csv_{sample}.log",
    params:
        seq_type="DNA",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    script:
        "../scripts/fasta_to_csv.py"


rule prots2csv:
    input:
        fa=rules.extract_prots.output.fa,
    output:
        csv=INT_SAMPLES_DIR / "annotation" / "{sample}" / "proteins.csv",
    log:
        LOGS / "samples" / "annotation" / "prots2csv_{sample}.log",
    params:
        seq_type="PROTEIN",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    script:
        "../scripts/fasta_to_csv.py"
