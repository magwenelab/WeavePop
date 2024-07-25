
# =================================================================================================
#   Per lineage | Extract CDS and protein sequences from reference genomes
# =================================================================================================
# Extract cds and protein sequences from reference genomes
rule extract_ref_seqs:
    input:
        gff = REFDIR / "{lineage}" / "{lineage}.gff",
        fasta = REFDIR / "{lineage}" / "{lineage}.fasta",
        config = rules.agat_config.output
    output:
        cds = REFDIR / "{lineage}" / "{lineage}.cds.fa",
        prots = REFDIR / "{lineage}" / "{lineage}.prots.fa"
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
        gff = REFDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "genes.gff",
        fasta = REFDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "sequences.fa",
        cds = REFDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "cds.fa",
        prots = REFDIR / "snpeff_data" / str(config["species_name"] + "_{lineage}") / "protein.fa"
    conda:
        "../envs/variants.yaml"
    params:
        name = config["species_name"] + "_{lineage}",
        config = REFDIR / "snpeff_data" / "snpEff.config",
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
        touch(REFDIR / "snpeff_data" / "{lineage}.done")
    conda:
        "../envs/variants.yaml"
    params:
        config = REFDIR/ "snpeff_data" / "snpEff.config",
        dir = os.getcwd() / REFDIR / "snpeff_data",
        name = config["species_name"] + "_{lineage}"
    log:
        "logs/references/build_dbs_{lineage}.log"
    shell:
        """
        snpEff build -gff3 -v -dataDir {params.dir} -config {params.config} {params.name} &>> {log}
        """
