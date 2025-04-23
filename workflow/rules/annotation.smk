# =================================================================================================
# Per sample | Run Liftoff to annotate the assembly with the corresponding reference genome
# =================================================================================================


rule liftoff:
    input:
        unpack(liftoff_input),
    output:
        ref_gff=INT_SAMPLES_DIR / "annotation"/ "liftoff" / "{sample}" / "ref.gff",
        gff=INT_SAMPLES_DIR / "annotation"/ "liftoff" / "{sample}" / "lifted.gff",
        polished=INT_SAMPLES_DIR / "annotation"/ "liftoff" / "{sample}" / "lifted.gff_polished",
        unmapped=INT_SAMPLES_DIR / "annotation"/ "liftoff" / "{sample}" / "unmapped_features.txt",
        intermediate=directory(INT_SAMPLES_DIR / "annotation"/ "liftoff" / "{sample}" / "intermediate_liftoff"),
        fai=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa.fai",
        mmi=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa.mmi",
    params:
        extra=config["liftoff"]["extra"],
        outpath=INT_SAMPLES_DIR / "annotation"/ "liftoff" / "{sample}",
    log:
        LOGS / "samples" / "annotation" / "liftoff_{sample}.log",
    threads: config["liftoff"]["threads"]
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
        gff=INT_SAMPLES_DIR / "annotation" / "{sample}" / "intergenic.gff",
    params:
        extra=config["agat"]["extra"],
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
        "{params.extra} "
        "&> {log} "


rule add_introns:
    input:
        gff=rules.add_intergenic.output.gff,
        config=rules.agat_config.output,
    output:
        gff=INT_SAMPLES_DIR / "annotation" / "{sample}" / "interg_introns.gff",
    params:
        extra=config["agat"]["extra"],
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
        "{params.extra} "
        "&> {log} "


rule annotation_gff2tsv:
    input:
        gff=rules.add_introns.output.gff,
        config=rules.agat_config.output,
    output:
        tsv=INT_SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff.tsv",
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
        fa=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa",
        config=rules.agat_config.output,
    output:
        fa=SAMPLES_DIR / "annotation" / "{sample}" / "cds.fa",
    params:
        extra=config["agat"]["extra"],
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
        "{params.extra} "
        "&> {log} "


rule extract_prots:
    input:
        gff=rules.reformat_annotation.output.gff,
        fa=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa",
        config=rules.agat_config.output,
        cds=rules.extract_cds.output.fa,
    output:
        fa=SAMPLES_DIR / "annotation" / "{sample}" / "proteins.fa",
    params:
        extra=config["agat"]["extra"],
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
        "-c {input.config}"
        "{params.extra} &> {log} "
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
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.fa} "
        "-s {wildcards.sample} "
        "-t DNA "
        "-o {output.csv} "
        "&> {log}"


rule prots2csv:
    input:
        fa=rules.extract_prots.output.fa,
    output:
        csv=INT_SAMPLES_DIR / "annotation" / "{sample}" / "proteins.csv",
    log:
        LOGS / "samples" / "annotation" / "prots2csv_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.fa} "
        "-s {wildcards.sample} "
        "-t PROTEIN "
        "-o {output.csv} "
        "&> {log}"
