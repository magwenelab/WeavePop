

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

rule cat_fastas:
    input:
        cds = expand(OUTDIR / "agat" / "{sample}/cds.fa",sample=SAMPLES),
        prots = expand(OUTDIR / "agat" / "{sample}/proteins.fa",sample=SAMPLES)
    output:
        cds = DATASET_OUTDIR / "cds.fa",
        prots = DATASET_OUTDIR / "proteins.fa"
    shell:
        "cat {input.cds} > {output.cds} "
        "&& "
        "cat {input.prots} > {output.prots}"

rule by_gene:
    input:
        cds = rules.cat_fastas.output.cds,
        prots = rules.cat_fastas.output.prots,
        ids = "protein_list.txt"
    output:
        cds = DATASET_OUTDIR / "cds.done",
        prots = DATASET_OUTDIR / "prots.done"
    params:
        cds = DATASET_OUTDIR / "by_cds",
        prots = DATASET_OUTDIR / "by_protein"
    script:
        "../scripts/by_cds.sh"