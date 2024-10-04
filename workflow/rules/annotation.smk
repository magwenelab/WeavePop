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


# rule move_liftoff:
#     input:
#         rules.liftoff.output.polished,
#     output:
#         SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff",
#     log:
#         "logs/samples/annotation/move_liftoff_{sample}.log",
#     conda:
#         "../envs/shell.yaml"
#     shell:
#         "mv {input} {output} 2> {log}"

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
        "logs/samples/annotation/add_intergenic_{sample}.log",
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
        gff=INT_SAMPLES_DIR / "annotation" / "{sample}" / "intergenic_and_introns.gff",
    params:
        extra=config["agat"]["extra"],
    log:
        "logs/samples/annotation/add_introns_{sample}.log",
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
        "logs/samples/annotation/annotation_gff2tsv_{sample}.log",
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

rule fix_annotation:
    input:
        tsv=rules.annotation_gff2tsv.output.tsv,
    output:
        tsv=SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff.tsv",
        gff=INT_SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff",
    log:
        "logs/samples/annotation/fix_annotation_{sample}.log",
    conda:
        "../envs/shell.yaml"
    script:
        "../scripts/fix_annotation.py"    

rule sort_gff:
    input:
        gff=rules.fix_annotation.output.gff,
        config=rules.agat_config.output,
    output:
        gff=SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff",
    log:
        "logs/samples/annotation/sort_{sample}.log",
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gxf2gxf.pl "
        "-g {input.gff} "
        "-o {output.gff} "
        "-c {input.config} "
        "&> {log} "


# =================================================================================================
# Per sample | Run AGAT to extract CDS, protein , intergenic and intronic sequences
# =================================================================================================


rule extract_cds:
    input:
        gff=rules.sort_gff.output.gff,
        fa=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa",
        config=rules.agat_config.output,
    output:
        fa=SAMPLES_DIR / "annotation" / "{sample}" / "cds.fa",
    params:
        extra=config["agat"]["extra"],
    log:
        "logs/samples/annotation/extract_cds_{sample}.log",
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
        gff=rules.sort_gff.output.gff,
        fa=SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa",
        config=rules.agat_config.output,
        cds=rules.extract_cds.output.fa,
    output:
        fa=SAMPLES_DIR / "annotation" / "{sample}" / "proteins.fa",
    params:
        extra=config["agat"]["extra"],
    log:
        "logs/samples/annotation/extract_prots_{sample}.log",
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
# Per sample | Convert fasta to csv
# =================================================================================================


rule cds2csv:
    input:
        fa=rules.extract_cds.output.fa,
    output:
        csv=INT_SAMPLES_DIR / "annotation" / "{sample}" / "cds.csv",
    log:
        "logs/samples/annotation/cds2csv_{sample}.log",
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
        "logs/samples/annotation/prots2db_{sample}.log",
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

