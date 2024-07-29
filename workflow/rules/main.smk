# =================================================================================================
# Per sample | Run Snippy to map reads to reference genome, get assembly and call SNPs
# =================================================================================================

def snippy_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "fq1": FQ_DATA / (s["sample"] + FQ1),
        "fq2": FQ_DATA / (s["sample"] + FQ2),
        "refgenome": s["refgenome"],
    }

rule snippy:
    input:
        unpack(snippy_input)
    output:
        fa = OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa",
        bam = OUTDIR / "snippy" / "{sample}" / "snps.bam",
        ref = OUTDIR / "snippy" / "{sample}" / "ref.fa",
        bai = OUTDIR / "snippy" / "{sample}" / "snps.bam.bai",
        vcf = OUTDIR / "snippy" / "{sample}" / "snps.vcf.gz"
    threads: 
        config["snippy"]["threads"]
    params:
        outpath = OUTDIR / "snippy",
        extra = config["snippy"]["extra"]
    conda:
        "../envs/snippy.yaml"
    log:
        "logs/samples/snippy/snippy_{sample}.log"
    shell:
        "snippy --outdir {params.outpath}/{wildcards.sample} "
        "--cpus {threads} "
        "--ref {input.refgenome} "
        "--R1 {input.fq1} "
        "--R2 {input.fq2} "
        "--force "
        "{params.extra} &> {log}"

# =================================================================================================
# Per sample | Run Liftoff to annotate the assembly with the coresponding reference genome
# =================================================================================================

def liftoff_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "target": OUTDIR / "snippy" / s["sample"] / "snps.consensus.fa" ,
        "refgff": s["refgff"],
        "refgenome": s["refgenome"],
    }
rule liftoff:
    input:
        unpack(liftoff_input)
    output:
        ref_gff = OUTDIR / "liftoff" / "{sample}" / "ref.gff",
        gff = OUTDIR / "liftoff" / "{sample}" / "lifted.gff",
        polished = OUTDIR / "liftoff" / "{sample}" / "lifted.gff_polished",
        unmapped = OUTDIR / "liftoff" / "{sample}" / "unmapped_features.txt",
        intermediate = directory(OUTDIR / "liftoff" / "{sample}" / "intermediate_files"),
        fai = OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa.fai",
        mmi = OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa.mmi"
    threads: 
        config["liftoff"]["threads"]
    params:
        extra = config["liftoff"]["extra"],
        outpath =  OUTDIR / "liftoff" / "{sample}"
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

# =================================================================================================
# Per sample | Run AGAT to extract CDS and protein sequences
# =================================================================================================

rule agat_cds:
    input:
        gff = rules.liftoff.output.polished,
        fa = rules.snippy.output.fa,
        config = rules.agat_config.output
    output:
        cds = OUTDIR / "agat" / "{sample}" / "cds.fa",
    conda:
        "../envs/agat.yaml"
    params:
        extra = config["agat"]["extra"]
    log: 
        "logs/samples/agat/agat_cds_{sample}.log",
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
        gff = rules.liftoff.output.polished,
        fa = rules.snippy.output.fa,
        config = rules.agat_config.output,
        cds = rules.agat_cds.output.cds
    output:
        prots = OUTDIR / "agat" / "{sample}" / "proteins.fa"
    conda:
        "../envs/agat.yaml"
    params:
        extra = config["agat"]["extra"]
    log: 
        "logs/samples/agat/agat_prots_{sample}.log"
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
        cds = OUTDIR / "agat" / "{sample}" / "cds.fa"
    output:
        cds = OUTDIR / "agat" / "{sample}" / "cds.csv"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/sequences/cds2csv_{sample}.log"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.cds} "
        "-s {wildcards.sample} "
        "-t DNA "
        "-o {output.cds} "
        "&> {log}"

rule prots2csv:
    input: 
        prots = OUTDIR / "agat" / "{sample}" / "proteins.fa"
    output:
        prots = OUTDIR / "agat" / "{sample}" / "proteins.csv"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/sequences/prots2db_{sample}.log"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.prots} "
        "-s {wildcards.sample} "
        "-t PROTEIN "
        "-o {output.prots} "
        "&> {log}"