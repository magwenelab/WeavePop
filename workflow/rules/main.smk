
rule snippy:
    input:
        unpack(snippy_input)
    output:
        fa = OUTDIR / "snippy" / "{sample}/snps.consensus.fa",
        bam = OUTDIR / "snippy" / "{sample}/snps.bam"
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
        "{params.extra} "
        "--force &> {log}"

rule liftoff:
    input:
        unpack(liftoff_input),
        features = FEATURE_FILE
    output:
        gff = OUTDIR / "liftoff" / "{sample}/lifted.gff",
        polished = OUTDIR / "liftoff" / "{sample}/lifted.gff_polished",
        unmapped = OUTDIR / "liftoff" / "{sample}/unmapped_features.txt" 
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
        "ln -s -r {input.refgff} {params.outpath}/{wildcards.sample}/ref.gff "
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
        "{input.refgenome} &> {log}"

rule agat:
    input:
        gff = rules.liftoff.output.gff,
        fa = rules.snippy.output.fa
    output:
        cds = OUTDIR / "agat" / "{sample}/predicted_cds.fa",
        prots = OUTDIR / "agat" / "{sample}/predicted_proteins.fa"
    conda:
        "../envs/agat.yaml"
    params:
        extra = config["agat"]["extra"]
    log: 
        cds = "logs/agat/{sample}_cds.log",
        prots = "logs/agat/{sample}_prots.log"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.cds} "
        "{params.extra} "
        "&> {log.cds} "
        " && "
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} " 
        "-f {input.fa} "
        "-o {output.prots} "
        "-p  &> {log.prots} "
        "{params.extra} " 
        " && "
        "rm lifted.agat.log || true"

rule agat_header:
    input:
        cds = rules.agat.output.cds,
        prots = rules.agat.output.prots
    output:
        cds = OUTDIR / "agat" / "{sample}/cds.fa",
        prots = OUTDIR / "agat" / "{sample}/proteins.fa"
    shell:
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' {input.cds} > {output.cds} "
        "&& "
        "seqkit replace -p '($)' -r ' sample={wildcards.sample}' {input.prots} > {output.prots} "
        "&& "
        "rm {input.cds} {input.prots}"