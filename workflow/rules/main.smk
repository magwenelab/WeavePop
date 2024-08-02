
# =================================================================================================
#   Tables for sample-dependent input files 
# =================================================================================================

def get_filtered_sample_reference():
    filtered_output = rules.quality_filter.output.metadata
    filtered_table = pd.read_csv(filtered_output, header=0)
    d={'sample': filtered_table["sample"],
        'lineage': filtered_table["lineage"],
        'refgenome': REFDIR / filtered_table["lineage"] / (filtered_table["lineage"] + ".fasta"),
        'refgff': REFDIR / filtered_table["lineage"] / (filtered_table["lineage"] + ".gff")}
    SAMPLE_REFERENCE = pd.DataFrame(data=d).set_index("sample", drop=False)
    LINEAGE_REFERENCE = pd.DataFrame(data=d).set_index("lineage", drop=False)
    return SAMPLE_REFERENCE, LINEAGE_REFERENCE

SAMPLE_REFERENCE, LINEAGE_REFERENCE = get_filtered_sample_reference()

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
        fa = OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa",
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
        fa = OUTDIR / "snippy" / "{sample}" / "snps.consensus.fa",
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