# Map reads to corresponding reference genome
rule snippy:
    input:
        unpack(snippy_input)
    output:
        fa = OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa",
        bam = OUTDIR / "snippy" / "{sample}" / "snps.bam",
        ref = OUTDIR / "snippy" / "{sample}" / "ref.fa"
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
rule agat_cds:
    input:
        gff = rules.liftoff.output.gff,
        fa = rules.snippy.output.fa
    output:
        cds = OUTDIR / "agat" / "{sample}" / "cds.fa",
    conda:
        "../envs/agat.yaml"
    params:
        extra = config["agat"]["extra"]
    log: 
        cds = "logs/agat/{sample}_cds.log",
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.cds} "
        "{params.extra} "
        "&> {log.cds} "

# Extract the amino acid sequence of each isoform of each gene of a sample
rule agat_prots:
    input:
        gff = rules.liftoff.output.gff,
        fa = rules.snippy.output.fa
    output:
        prots = OUTDIR / "agat" / "{sample}" / "proteins.fa"
    conda:
        "../envs/agat.yaml"
    params:
        extra = config["agat"]["extra"]
    log: 
        prots = "logs/agat/{sample}_prots.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.prots} "
        "-p  &> {log.prots} "
        "{params.extra} " 
        "&& "
        "sed -i 's/type=cds//g' {output} &>> {log.prots} "
        "&& "
        "rm lifted.agat.log &>> {log.prots} || true"

rule by_id:
    input: 
        prots = rules.agat_prots.output.prots,
        cds = rules.agat_cds.output.cds,
    output:
        touch(OUTDIR / "agat" / "{sample}" / "cds.done"),
        touch(OUTDIR / "agat" / "{sample}" / "prots.done")
    params:
        db = DATASET_OUTDIR / "sequences.db"
    log:
        "logs/agat/database_{sample}.log"
    shell:
        "python workflow/scripts/sqlfasta.py populate-db {params.db} {input.prots} {wildcards.sample} -t Protein &> {log}"
        "&& "
        "python workflow/scripts/sqlfasta.py populate-db {params.db} {input.cds} {wildcards.sample} -t DNA &>> {log}"
