# =================================================================================================
#   Per lineage | Standardize GFF format and convert to TSV
# =================================================================================================

# Run AGAT to add and modify tags and convert to TSV 
rule agat_fix_gff:
    input:
        gff = REF_DATA / "{lineage}.gff",
        config = rules.agat_config.output
    output:
        fixed_ID = temp(REFS_DIR / "{lineage}" / "{lineage}.fixed.gff"),
        fixed_locus = temp(REFS_DIR / "{lineage}" / "{lineage}.fixed_locus.gff"),
        fixed_description = temp(REFS_DIR / "{lineage}" / "{lineage}.fixed_description.gff"),
        tsv = temp(REFS_DIR / "{lineage}" / "{lineage}.tsv")
    log:
        "logs/references/agat_fix_gff_{lineage}.log"
    resources:
        tmpdir = TEMPDIR
    conda:
        "../envs/agat.yaml"
    shell:
        """
        agat_convert_sp_gxf2gxf.pl -g {input.gff} -o {output.fixed_ID} -c {input.config} &> {log} && 
        agat_sq_add_locus_tag.pl --gff {output.fixed_ID} --li ID -o {output.fixed_locus} -c {input.config} &>> {log} && 
        agat_sp_manage_attributes.pl --gff {output.fixed_locus} --tag product/description -o {output.fixed_description} -c {input.config} &>> {log} && 
        agat_convert_sp_gff2tsv.pl --gff {output.fixed_description} -o {output.tsv} -c {input.config} &>> {log} 
        """

# Recreate IDs
rule fix_gff_tsv:
    input:
        tsv = rules.agat_fix_gff.output.tsv
    output:
        gff = REFS_DIR / "{lineage}.gff",
        tsv = REFS_DIR / "{lineage}" / "{lineage}.gff.tsv"
    log:
        "logs/references/fix_gff_tsv_{lineage}.log"
    resources:
        tmpdir = TEMPDIR
    shell:
        """
        python workflow/scripts/fix_gff.py -i {input.tsv} -og {output.gff} -ot {output.tsv} &> {log}
        """

