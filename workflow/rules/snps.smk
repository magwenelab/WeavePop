
# =================================================================================================
#   Per lineage | Extract CDS and protein sequences from reference genomes
# =================================================================================================
# Extract cds and protein sequences from reference genomes
rule extract_ref_seqs:
    input:
        gff = REFS_DIR / "{lineage}.gff",
        fasta = INT_REFS_DIR / "{lineage}" / "{lineage}.fasta",
        config = rules.agat_config.output
    output:
        cds = INT_REFS_DIR / "{lineage}" / "{lineage}.cds.fa",
        prots = INT_REFS_DIR / "{lineage}" / "{lineage}.prots.fa"
    conda:
        "../envs/agat.yaml"
    log:
        "logs/references/extract_ref_seqs_{lineage}.log"
    shell:
        "agat_sp_extract_sequences.pl -g {input.gff} -f {input.fasta} -o {output.cds} -c {input.config} &> {log} "
        "&& "
        "agat_sp_extract_sequences.pl -g {input.gff} -f {input.fasta} -o {output.prots} -p -c {input.config} &>> {log} "

# Make symbolic links in the snpeff_data directory and create config file
rule prepare_refs_db:
    input:
        gff = rules.extract_ref_seqs.input.gff,
        fasta = rules.extract_ref_seqs.input.fasta,
        cds = rules.extract_ref_seqs.output.cds,
        prots = rules.extract_ref_seqs.output.prots
    output:
        gff = INT_REFS_DIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "genes.gff",
        fasta = INT_REFS_DIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "sequences.fa",
        cds = INT_REFS_DIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "cds.fa",
        prots = INT_REFS_DIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "protein.fa"
    conda:
        "../envs/variants.yaml"
    params:
        name = config["species_name"] + "_{lineage}",
        config = INT_REFS_DIR / "snpeff_data" / "snpEff.config",
    log:
        "logs/references/prepare_dbs_{lineage}.log"
    shell:
        """
        echo "{params.name}.genome : {params.name}" >> {params.config} 2> {log} && 
        ln -s -r {input.gff} {output.gff} &>> {log} && 
        ln -s -r {input.fasta} {output.fasta} &>> {log} && 
        ln -s -r {input.cds} {output.cds} &>> {log} && 
        ln -s -r {input.prots} {output.prots} &>> {log} 
        """

# Build snpeff database for the reference genomes
rule build_refs_db:
    input:
        gff = rules.prepare_refs_db.output.gff,
        fasta = rules.prepare_refs_db.output.fasta,
        cds = rules.prepare_refs_db.output.cds,
        prots = rules.prepare_refs_db.output.prots
    output:
        touch(INT_REFS_DIR / "snpeff_data" / "{lineage}.done")
    conda:
        "../envs/variants.yaml"
    params:
        config = INT_REFS_DIR/ "snpeff_data" / "snpEff.config",
        dir = os.getcwd() / INT_REFS_DIR / "snpeff_data",
        name = config["species_name"] + "_{lineage}"
    log:
        "logs/references/build_dbs_{lineage}.log"
    shell:
        """
        snpEff build -gff3 -v -dataDir {params.dir} -config {params.config} {params.name} &>> {log}
        """

# =================================================================================================
#   Join samples per lineage | Intersect VCF files of all samples from each lineage and annotate
# =================================================================================================

def intersect_vcfs_input(wildcards):
    sample_wildcards = listing_samples(wildcards)
    l = LINEAGE_REFERENCE[LINEAGE_REFERENCE["sample"].isin(sample_wildcards)] 
    l = l.loc[wildcards.lineage,]
    return {
        "vcfs" : expand(SAMPLES_DIR / "snippy" / "{sample}" / "snps.vcf.gz", sample=l["sample"])
    }
    
rule intersect_vcfs:
    input:
        unpack(intersect_vcfs_input)
    output:
        vcf = INT_DATASET_DIR / "snps" / "{lineage}_intersection.vcf",
        tsv = INT_DATASET_DIR / "snps" / "{lineage}_presence.tsv"
    params:
        tmp_dir = "tmp_{lineage}"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/intersect_vcfs_{lineage}.log"
    shell:
        "xonsh workflow/scripts/intersect_vcfs.xsh "
        "-v {output.vcf} "
        "-p {output.tsv} "
        "-l {wildcards.lineage} "
        "-t {params.tmp_dir} "
        "{input.vcfs} "
        "&> {log}"

rule snpeff:
    input:
        vcf = rules.intersect_vcfs.output.vcf,
        db_done = rules.build_refs_db.output
    output:
        vcf = INT_DATASET_DIR / "snps" / "{lineage}_snpeff.vcf",
        html = INT_DATASET_DIR / "snps" / "{lineage}_snpeff.html"
    params:
       dir = os.getcwd() / INT_REFS_DIR / "snpeff_data",
       config = INT_REFS_DIR / "snpeff_data" / "snpEff.config",
       name = config["species_name"] + "_{lineage}"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/snpeff_{lineage}.log"
    shell:
        "snpEff ann -v -classic "
        "-dataDir {params.dir} "
        "-config {params.config} "
        "-s {output.html} "
        "{params.name} "
        "{input.vcf} "
        "1> {output.vcf} 2> {log}"

rule extract_vcf_annotation:
    input:
        vcf = rules.snpeff.output.vcf
    output:
        effects = INT_DATASET_DIR / "snps" / "{lineage}_effects.tsv",
        variants = INT_DATASET_DIR / "snps" / "{lineage}_variants.tsv",
        lofs = INT_DATASET_DIR / "snps" / "{lineage}_lofs.tsv",
        nmds = INT_DATASET_DIR / "snps" / "{lineage}_nmds.tsv"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/dataset/snps/extract_vcf_annotation_{lineage}.log"
    shell:
        "xonsh workflow/scripts/extract_vcf_annotation.xsh "
        "-i {input.vcf} "
        "-e {output.effects} "
        "-v {output.variants} "
        "-f {output.lofs} "
        "-n {output.nmds} "
        "-l {wildcards.lineage} "
        "&> {log}"

