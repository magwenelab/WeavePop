# =================================================================================================
# Per sample | Run Liftoff to annotate the assembly with the coresponding reference genome
# =================================================================================================


rule liftoff:
    input:
        unpack(liftoff_input),
    output:
        ref_gff=INT_SAMPLES_DIR / "liftoff" / "{sample}" / "ref.gff",
        gff=INT_SAMPLES_DIR / "liftoff" / "{sample}" / "lifted.gff",
        polished=INT_SAMPLES_DIR / "liftoff" / "{sample}" / "lifted.gff_polished",
        unmapped=INT_SAMPLES_DIR / "liftoff" / "{sample}" / "unmapped_features.txt",
        intermediate=directory(INT_SAMPLES_DIR / "liftoff" / "{sample}" / "intermediate_liftoff"),
        fai=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa.fai",
        mmi=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa.mmi",
    params:
        extra=config["liftoff"]["extra"],
        outpath=INT_SAMPLES_DIR / "liftoff" / "{sample}",
    log:
        "logs/samples/liftoff/liftoff_{sample}.log",
    threads: config["liftoff"]["threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/liftoff.yaml"
    message:
        "Running Liftoff for {wildcards.sample}"
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


rule move_liftoff:
    input:
        rules.liftoff.output.polished,
    output:
        SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff",
    log:
        "logs/samples/annotation/move_liftoff_{sample}.log",
    conda:
        "../envs/shell.yaml"
    shell:
        "mv {input} {output} 2> {log}"


# =================================================================================================
# Per sample | Run AGAT to extract CDS and protein sequences
# =================================================================================================


rule agat_cds:
    input:
        gff=rules.move_liftoff.output,
        fa=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa",
        config=rules.agat_config.output,
    output:
        cds=SAMPLES_DIR / "annotation" / "{sample}" / "cds.fa",
    params:
        extra=config["agat"]["extra"],
    log:
        "logs/samples/annotation/agat_cds_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} "
        "-f {input.fa} "
        "-o {output.cds} "
        "-c {input.config} "
        "{params.extra} "
        "&> {log} "


rule agat_prots:
    input:
        gff=rules.move_liftoff.output,
        fa=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa",
        config=rules.agat_config.output,
        cds=rules.agat_cds.output.cds,
    output:
        prots=SAMPLES_DIR / "annotation" / "{sample}" / "proteins.fa",
    params:
        extra=config["agat"]["extra"],
    log:
        "logs/samples/annotation/agat_prots_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} "
        "-f {input.fa} "
        "-o {output.prots} "
        "-p "
        "-c {input.config}"
        "{params.extra} &> {log} "
        "&& "
        "sed -i 's/type=cds//g' {output} &>> {log} "


# =================================================================================================
# Per sample | Convert fasta to csv
# =================================================================================================


rule cds2csv:
    input:
        cds=SAMPLES_DIR / "annotation" / "{sample}" / "cds.fa",
    output:
        cds=INT_SAMPLES_DIR / "annotation" / "{sample}" / "cds.csv",
    log:
        "logs/samples/annotation/cds2csv_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.cds} "
        "-s {wildcards.sample} "
        "-t DNA "
        "-o {output.cds} "
        "&> {log}"


rule prots2csv:
    input:
        prots=SAMPLES_DIR / "annotation" / "{sample}" / "proteins.fa",
    output:
        prots=INT_SAMPLES_DIR / "annotation" / "{sample}" / "proteins.csv",
    log:
        "logs/samples/annotation/prots2db_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.prots} "
        "-s {wildcards.sample} "
        "-t PROTEIN "
        "-o {output.prots} "
        "&> {log}"

# =================================================================================================
# Per sample | Annotate and extract intergenic regions
# =================================================================================================


rule agat_intergenic:
    input:
        gff=rules.move_liftoff.output,
        config=rules.agat_config.output,
    output:
        gff=SAMPLES_DIR / "annotation" / "{sample}" / "intergenic.gff",
    params:
        extra=config["agat"]["extra"],
    log:
        "logs/samples/annotation/agat_intergenic_{sample}.log",
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

rule agat_extract_intergenic:
    input:
        gff=rules.agat_intergenic.output.gff,
        fa=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa",
        config=rules.agat_config.output,
    output:
        fa=SAMPLES_DIR / "annotation" / "{sample}" / "intergenic.fa",
    params:
        extra=config["agat"]["extra"],
    log:
        "logs/samples/annotation/agat_intergenic_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} "
        "-f {input.fa} "
        "-t intergenic_region "
        "-o {output.fa} "
        "-c {input.config} "
        "{params.extra} "
        "&> {log} "

rule intergenic2csv:
    input:
        inter=SAMPLES_DIR / "annotation" / "{sample}" / "intergenic.fa",
    output:
        inter=INT_SAMPLES_DIR / "annotation" / "{sample}" / "intergenic.csv",
    log:
        "logs/samples/annotation/intergenic2csv_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.inter} "
        "-s {wildcards.sample} "
        "-t INTERGENIC "
        "-o {output.inter} "
        "&> {log}"

rule intergenic_gff2tsv:
    input:
        gff=rules.agat_intergenic.output.gff,
        config=rules.agat_config.output,
    output:
        tsv=INT_SAMPLES_DIR / "annotation" / "{sample}" / "intergenic.gff.tsv",
    log:
        "logs/samples/annotation/intergenic_gff2tsv_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "--gff {input.gff} "
        "-o {output.tsv} "
        "-c {input.config} "
        "&> {log} "

rule intergenic_coords:
    input:
        csv=rules.intergenic2csv.output,
        gff=rules.intergenic_gff2tsv.output,
    output:
        csv=INT_SAMPLES_DIR / "annotation" / "{sample}" / "intergenic_coords.csv",
    log:
        "logs/samples/annotation/intergenic_coords_{sample}.log",
    conda:
        "../envs/shell.yaml"
    script:
        "../scripts/intergenic_coords.py"