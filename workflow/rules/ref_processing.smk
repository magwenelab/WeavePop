# ==================================================================================================
#   Per lineage | Add repetitive sequences, introns, intergenic regions, and convert to TSV
# ==================================================================================================
rule ref_add_intergenic:
    input:
        gff=INT_REFS_DIR / "{lineage}" / "{lineage}_annotated.gff",
        config=rules.agat_config.output,
    output:
        INT_REFS_DIR / "{lineage}" / "{lineage}_intergenic.gff",
    params:
        extra=config["agat"]["extra"],
    log:
        LOGS / "references" / "annotation" / "ref_add_intergenic_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_add_intergenic_regions.pl "
        "-g {input.gff} "
        "-o {output} "
        "-c {input.config} "
        "{params.extra} "
        "&> {log}"

rule ref_add_introns:
    input:
        gff=rules.ref_add_intergenic.output,
        config=rules.agat_config.output,
    output:
        gff=INT_REFS_DIR / "{lineage}" / "{lineage}_interg_introns.gff",
    params:
        extra=config["agat"]["extra"],
    log:
        LOGS / "references" / "ref_processing" / "ref_add_introns_{lineage}.log",
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
        "&> {log}"


rule ref_gff2tsv:
    input:
        target=rules.ref_add_introns.output.gff,
        config=rules.agat_config.output,
    output:
        tsv=INT_REFS_DIR / "{lineage}" / "{lineage}_interg_introns.gff.tsv",
    log:
        LOGS / "references" / "ref_processing" / "gff2tsv_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "-gff {input.target} "
        "-c {input.config} "
        "-o {output.tsv} "
        "&> {log}"

rule ref_add_repeats:
    input:
        gff=rules.ref_add_introns.output.gff,
        repeats=REFS_DIR / "{lineage}" / "{lineage}_repeats.bed",
    output:
        INT_REFS_DIR / "{lineage}" / "{lineage}_repeats.gff",
    log:
        LOGS / "references" / "ref_processing" / "ref_add_repeats_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/ref_add_repeats.xsh "
        "-g {input.gff} "
        "-r {input.repeats} "
        "-o {output} "
        "&> {log}"


rule ref_gff2tsv_2:
    input:
        target=rules.ref_add_repeats.output,
        config=rules.agat_config.output,
    output:
        tsv=INT_REFS_DIR / "{lineage}" / "{lineage}_repeats.gff.tsv",
    log:
        LOGS / "references" / "ref_processing" / "gff2tsv_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_convert_sp_gff2tsv.pl "
        "-gff {input.target} "
        "-c {input.config} "
        "-o {output.tsv} "
        "&> {log}"


rule ref_reformat_annotation:
    input:
        tsv=rules.ref_gff2tsv_2.output.tsv,
    output:
        tsv=REFS_DIR / "{lineage}" / "{lineage}.gff.tsv",
        gff=REFS_DIR / "{lineage}" / "{lineage}.gff",
    params:
        lineage="{lineage}",
        version="lineage",
    log:
        LOGS / "references" / "ref_processing" / "ref_reformat_annotation_{lineage}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/reformat_annotation.py"

# =================================================================================================
#   Per lineage | Extract CDS and protein sequences from reference genomes and convert to CSV
# =================================================================================================


rule extract_cds_seqs:
    input:
        gff=INT_REFS_DIR / "{lineage}" / "{lineage}_interg_introns.gff",
        fasta=INT_REFS_DIR / "{lineage}" / "{lineage}.fasta",
        config=rules.agat_config.output,
    output:
        cds=INT_REFS_DIR / "{lineage}" / "{lineage}.cds.fa",
    log:
        LOGS / "references" / "snpeff" / "extract_cds_seqs_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} "
        "-f {input.fasta} "
        "-o {output.cds} "
        "-c {input.config} "
        "&> {log}"

rule extract_protein_seqs:
    input:
        gff=INT_REFS_DIR / "{lineage}" / "{lineage}_interg_introns.gff",
        fasta=INT_REFS_DIR / "{lineage}" / "{lineage}.fasta",
        config=rules.agat_config.output,
        cds=rules.extract_cds_seqs.output.cds,
    output:
        prots=INT_REFS_DIR / "{lineage}" / "{lineage}.prots.fa",
    log:
        LOGS / "references" / "snpeff" / "extract_protein_seqs_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl "
        "-g {input.gff} "
        "-f {input.fasta} "
        "-o {output.prots} "
        "-p "
        "-c {input.config} "
        "&> {log}"

rule ref_cds2csv:
    input:
        fa=INT_REFS_DIR / "{lineage}" / "{lineage}.cds.fa",
    output:
        csv=INT_REFS_DIR / "{lineage}" / "{lineage}.cds.csv",
    log:
        LOGS / "references" / "annotation" / "ref_cds2csv_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.fa} "
        "-l {wildcards.lineage} "
        "-t DNA "
        "-o {output.csv} "
        "&> {log}"

rule ref_prots2csv:
    input:
        fa=INT_REFS_DIR / "{lineage}" / "{lineage}.prots.fa",
    output:
        csv=INT_REFS_DIR / "{lineage}" / "{lineage}.prots.csv",
    log:
        LOGS / "references" / "annotation" / "ref_prots2csv_{lineage}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/variants.yaml"
    shell:
        "python workflow/scripts/fasta_to_csv.py "
        "-f {input.fa} "
        "-l {wildcards.lineage} "
        "-t PROTEIN "
        "-o {output.csv} "
        "&> {log}"

# =================================================================================================
#   All refernces | Obtain chromosome lengths
# =================================================================================================


rule chromosome_lengths:
    input:
        INT_REFS_DIR / "{lineage}" / "{lineage}.fasta",         
    output:
        INT_REFS_DIR / "{lineage}" / "chromosome_lengths.tsv",
    log:
        LOGS / "references" / "ref_processing" / "{lineage}_chromosome_lengths.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/agat.yaml"
    shell:
        """
        seqkit fx2tab -l -i -n {input} |\
        awk -v lin={wildcards.lineage} '{{print lin, $0}}' OFS='\t' \
        1> {output} 2> {log}
        """

rule join_chromosome_lengths:
    input:
        chrom_names=CHROM_NAMES,  
        chrom_lengths=expand(INT_REFS_DIR / "{lineage}" / "chromosome_lengths.tsv", lineage=LINEAGES),
    output:
        INT_REFS_DIR / "chromosome_lengths.tsv",
    log:
        LOGS / "references" / "ref_processing" / "join_chromosome_lengths.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_chromosome_lengths.py"        