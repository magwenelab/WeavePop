# =================================================================================================
# Per sample | Run Liftoff to annotate the assembly with the coresponding reference genome
# =================================================================================================

def liftoff_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "target": SAMPLES_DIR / "snippy" / s["sample"] / "snps.consensus.fa" ,
        "refgff": s["refgff"],
        "refgenome": s["refgenome"],
    }
rule liftoff:
    input:
        unpack(liftoff_input)
    output:
        ref_gff = INT_SAMPLES_DIR / "liftoff" / "{sample}" / "ref.gff",
        gff = INT_SAMPLES_DIR / "liftoff" / "{sample}" / "lifted.gff",
        polished = INT_SAMPLES_DIR / "liftoff" / "{sample}" / "lifted.gff_polished",
        unmapped = INT_SAMPLES_DIR / "liftoff" / "{sample}" / "unmapped_features.txt",
        intermediate = directory(INT_SAMPLES_DIR / "liftoff" / "{sample}" / "intermediate_liftoff"),
        fai = SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa.fai",
        mmi = SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa.mmi"
    threads: 
        config["liftoff"]["threads"]
    params:
        extra = config["liftoff"]["extra"],
        outpath =  INT_SAMPLES_DIR / "liftoff" / "{sample}"
    conda:
        "../envs/liftoff.yaml"
    log:
        "logs/samples/liftoff/liftoff_{sample}.log" 
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
        rules.liftoff.output.polished
    output:
        SAMPLES_DIR / "annotation" / "{sample}" / "annotation.gff"
    shell:
        "mv {input} {output}"
# =================================================================================================
# Per sample | Run AGAT to extract CDS and protein sequences
# =================================================================================================

rule agat_cds:
    input:
        gff = rules.move_liftoff.output,
        fa = SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa",
        config = rules.agat_config.output
    output:
        cds = SAMPLES_DIR / "annotation" / "{sample}" / "cds.fa",
    conda:
        "../envs/agat.yaml"
    params:
        extra = config["agat"]["extra"]
    log: 
        "logs/samples/annotation/agat_cds_{sample}.log",
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
        gff = rules.move_liftoff.output,
        fa = SAMPLES_DIR / "snippy" / "{sample}" / "snps.consensus.fa",
        config = rules.agat_config.output,
        cds = rules.agat_cds.output.cds
    output:
        prots = SAMPLES_DIR / "annotation" / "{sample}" / "proteins.fa"
    conda:
        "../envs/agat.yaml"
    params:
        extra = config["agat"]["extra"]
    log: 
        "logs/samples/annotation/agat_prots_{sample}.log"
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
        cds = SAMPLES_DIR / "annotation" / "{sample}" / "cds.fa"
    output:
        cds = INT_SAMPLES_DIR / "annotation" / "{sample}" / "cds.csv"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/samples/annotation/cds2csv_{sample}.log"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.cds} "
        "-s {wildcards.sample} "
        "-t DNA "
        "-o {output.cds} "
        "&> {log}"

rule prots2csv:
    input: 
        prots = SAMPLES_DIR / "annotation" / "{sample}" / "proteins.fa"
    output:
        prots = INT_SAMPLES_DIR / "annotation" / "{sample}" / "proteins.csv"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/samples/annotation/prots2db_{sample}.log"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.prots} "
        "-s {wildcards.sample} "
        "-t PROTEIN "
        "-o {output.prots} "
        "&> {log}"