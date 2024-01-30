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
        "{params.extra} "
        "--force &> {log}"

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
        " && "
        "rm lifted.agat.log || true"

# checkpoint mock:
#     output:
#         directory(DATASET_OUTDIR / "cds")
#     shell:
#         "mkdir -p {output} && "
#         "echo {params} && "
#         "cut -f12 results/references/FungiDB-65_CneoformansH99.tsv | sort | uniq | grep CNAG | while read line; do touch {output}/$line.fa; echo someline2 >> {output}/$line.fa ; done"

# Make a fasta file for each isoform with the nucleotide sequence of samples it is present in
rule by_id_cds:
    input:
        rules.agat_cds.output.cds
    output:
        done = touch(OUTDIR / "agat" / "{sample}" / "by_cds.done")
    params: 
        script = workflow.source_path("../scripts/by_id.py"), 
        outdir = directory(DATASET_OUTDIR / "cds")
    shell:
        "python {params.script} {input} {wildcards.sample} --outdir {params.outdir} "

# Make a fasta file for each isoform with the amino acid sequence of samples it is present in
rule by_id_proteins:
    input:
        rules.agat_prots.output.prots
    output:
        touch(OUTDIR / "agat" / "{sample}" / "by_proteins.done")
    params: 
        script = workflow.source_path("../scripts/by_id.py"),
        outdir = DATASET_OUTDIR / "proteins"                 
    shell:
        "python {params.script} {input} {wildcards.sample} --outdir {params.outdir} "