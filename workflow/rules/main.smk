# Map reads to corresponding reference genome
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
        "logs/snippy/{sample}.log"
    shell:
        "snippy --outdir {params.outpath}/{wildcards.sample} "
        "--cpus {threads} "
        "--ref {input.refgenome} "
        "--R1 {input.fq1} "
        "--R2 {input.fq2} "
        "--force "
        "{params.extra} &> {log}"

# Lift over annotation from corresponding reference genome to sample assembly
def liftoff_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "target": OUTDIR / "snippy" / s["sample"] / "snps.consensus.fa" ,
        "refgff": s["refgff"],
        "refgenome": s["refgenome"],
    }
rule liftoff:
    input:
        unpack(liftoff_input),
        features = FEATURE_FILE
    output:
        gff = OUTDIR / "liftoff" / "{sample}" / "lifted.gff",
        polished = OUTDIR / "liftoff" / "{sample}" / "lifted.gff_polished",
        unmapped = OUTDIR / "liftoff" / "{sample}" / "unmapped_features.txt" 
    threads: 
        config["liftoff"]["threads"]
    params:
        extra = config["liftoff"]["extra"],
        outpath =  OUTDIR / "liftoff",
    conda:
        "../envs/liftoff.yaml"
    log:
        "logs/liftoff/{sample}.log" 
    shell:
        "ln -s -r {input.refgff} {params.outpath}/{wildcards.sample}/ref.gff &> {log} || true "
        "&& "
        "liftoff "
        "-g {params.outpath}/{wildcards.sample}/ref.gff " 
        "-f {input.features} "
        "-o {output.gff} "
        "-dir {params.outpath}/{wildcards.sample}/intermediate_files "
        "-u {params.outpath}/{wildcards.sample}/unmapped_features.txt "
        "-p {threads} "
        "-polish "
        "{params.extra} "
        "{input.target} "
        "{input.refgenome} &>> {log}"

# Extract the nucleotide sequence of each isoform of each gene of a sample
rule agat_config:
    output:
        "config/agat_config.yaml"
    conda:
        "../envs/agat.yaml"
    log:
        "logs/agat_config.log"
    shell:
        "agat config --expose 2> {log} && "
        "mv agat_config.yaml {output} 2> {log} && "
        "sed -i 's/log: true/log: false/g' {output} &>> {log} "

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
        "logs/agat/cds_{sample}.log",
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.cds} "
        "-c {input.config} "
        "{params.extra} "
        "&> {log} "

# Extract the amino acid sequence of each isoform of each gene of a sample
rule agat_prots:
    input:
        gff = rules.liftoff.output.polished,
        fa = rules.snippy.output.fa,
        config = rules.agat_config.output
    output:
        prots = OUTDIR / "agat" / "{sample}" / "proteins.fa"
    conda:
        "../envs/agat.yaml"
    params:
        extra = config["agat"]["extra"]
    log: 
        "logs/agat/prots_{sample}.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.prots} "
        "-p "
        "-c config/agat_config.yaml "
        "{params.extra} &> {log} " 
        "&& "
        "sed -i 's/type=cds//g' {output} &>> {log} "

# Make SQL database with cds of all samples
rule cds2db:
    input: 
        cds = rules.agat_cds.output.cds,
    output:
        touch(DATASET_OUTDIR / "database" / "{sample}" / "cds.done"),
    params:
        db = DATASET_OUTDIR / "sequences.db"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/cds2db_{sample}.log"
    shell:
        "python workflow/scripts/build_sequences_db.py -d {params.db} -f {input.cds} -s {wildcards.sample} -t DNA &> {log}"

# Make SQL database with proteins of all samples
rule prots2db:
    input: 
        prots = rules.agat_prots.output.prots,
    output:
        touch(DATASET_OUTDIR / "database" / "{sample}" / "prots.done")
    params:
        db = DATASET_OUTDIR / "sequences.db"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/prots2db_{sample}.log"
    shell:
        "python workflow/scripts/build_sequences_db.py -d {params.db} -f {input.prots} -s {wildcards.sample} -t PROTEIN &> {log}"