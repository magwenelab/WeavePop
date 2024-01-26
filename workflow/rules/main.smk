
# rule samples_list:
#     output: 
#         "files/samples.txt"
#     run:
#         sample = samplefile["sample"]
#         sample.to_csv("files/samples.txt", index = False, header = False)

# rule ref_agat:
#     input: 
#         lin_gff = REFDIR + "{lineage}.gff",
#         lin_fasta = REFDIR + "{lineage}.fasta"
#     output:
#         cds = REFDIR + "{lineage}_predicted_cds.fa",
#         prots = REFDIR + "{lineage}_predicted_proteins.fa"
#     conda:
#         "../envs/agat.yaml"
#     log:
#         cds = "logs/references/{lineage}_ref_agat_cds.log",   
#         prots = "logs/references/{lineage}_ref_agat_prots.log"
#     shell:
#         "agat_sp_extract_sequences.pl "
#         "-g {input.lin_gff} "
#         "-f {input.lin_fasta} "
#         "-o {output.cds} "
#         "&> {log.cds} "
#         " && "
#         "agat_sp_extract_sequences.pl "
#         "-g {input.lin_gff} "
#         "-f {input.lin_fasta} "
#         "-o {output.prots} "
#         "-p &> {log.prots} "
#         " && "
#         "rm {wildcards.lineage}.agat.log || true"  

# rule protein_list:
#     input:
#         fasta = REFDIR + "{lineage}_predicted_proteins.fa"
#     output:
#         list = REFDIR + "{lineage}_protein_list.txt"
#     conda:
#         "../envs/agat.yaml"
#     log:
#         "logs/references/{lineage}_protein_list.log"
#     shell:
#         "seqkit seq -n -i {input.fasta} 1> {output.list} 2> {log}"

# rule cat_lists:
#     input: 
#         expand(REFDIR + "{lineage}_protein_list.txt", lineage=LINS)
#     output:
#         "files/protein_list.txt"
#     log:
#         "logs/references/cat_list.log"
#     shell:
#         "cat {input} | sort | uniq > {output} 2> {log}"

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

# rule agat:
#     input:
#         gff = "analysis/{sample}/lifted.gff_polished",
#         fa = "analysis/{sample}/snps.consensus.fa"
#     output:
#         cds = "analysis/{sample}/predicted_cds.fa",
#         prots = "analysis/{sample}/predicted_proteins.fa"
#     conda:
#         "../envs/agat.yaml"
#     log: 
#         cds = "logs/agat/{sample}_cds.log",
#         prots = "logs/agat/{sample}_prots.log"
#     shell:
#         "agat_sp_extract_sequences.pl "
#         "-g {input.gff} " 
#         "-f {input.fa} "
#         "-o {output.cds} "
#         "&> {log.cds} "
#         " && "
#         "agat_sp_extract_sequences.pl "
#         "-g {input.gff} " 
#         "-f {input.fa} "
#         "-o {output.prots} "
#         "-p  &> {log.prots} " 
#         " && "
#         "rm lifted.agat.log || true"

# rule index_proteins:
#     input:
#         "analysis/{sample}/predicted_proteins.fa"
#     output:
#         "analysis/{sample}/predicted_proteins.fa.fai"
#     conda:
#         "../envs/agat.yaml"        
#     log:
#         "logs/faidx/{sample}_proteins.log"    
#     shell:
#         "seqkit faidx {input} &> {log}"

# rule index_cds:
#     input:
#         "analysis/{sample}/predicted_cds.fa"
#     output:
#         "analysis/{sample}/predicted_cds.fa.fai"
#     conda:
#         "../envs/agat.yaml"        
#     log:
#         "logs/faidx/{sample}_cds.log"      
#     shell:
#         "seqkit faidx {input} &> {log}"

# rule by_proteins:
#     input:
#         "files/protein_list.txt",
#         "files/samples.txt",
#         expand("analysis/{sample}/predicted_proteins.fa",sample=SAMPLES),
#         expand("analysis/{sample}/predicted_proteins.fa.fai",sample=SAMPLES)
#     output:
#         "results/proteins.done"
#     log:
#         "logs/proteins/proteins.log"  
#     script:
#         "scripts/by_proteins.sh"

# rule by_cds:
#     input:
#         "files/protein_list.txt",
#         "files/samples.txt",
#         expand("analysis/{sample}/predicted_cds.fa",sample=SAMPLES),
#         expand("analysis/{sample}/predicted_cds.fa.fai",sample=SAMPLES)
#     output:
#         "results/cds.done"
#     log:
#         "logs/cds/cds.log"  
#     script:
#         "scripts/by_cds.sh"
