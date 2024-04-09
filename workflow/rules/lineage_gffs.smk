rule agat_fix_gff:
    input:
        gff = REF_DATA / "{lineage}.gff"
    output:
        fixed_ID = temp(REFDIR / "{lineage}" / "{lineage}.fixed.gff"),
        fixed_locus = temp(REFDIR / "{lineage}" / "{lineage}.fixed_locus.gff"),
        fixed_description = temp(REFDIR / "{lineage}" / "{lineage}.fixed_description.gff"),
        tsv = temp(REFDIR / "{lineage}" / "{lineage}.tsv")
    conda:
        "../envs/agat.yaml"
    log:
        "logs/references/agat_fix_gff_{lineage}.log"
    shell:
        """
        agat_convert_sp_gxf2gxf.pl -g {input.gff} -o {output.fixed_ID} &> {log} && 
        agat_sq_add_locus_tag.pl --gff {output.fixed_ID} --li ID -o {output.fixed_locus} &>> {log} && 
        agat_sp_manage_attributes.pl --gff {output.fixed_locus} --tag product/description -o {output.fixed_description} &>> {log} && 
        agat_convert_sp_gff2tsv.pl --gff {output.fixed_description} -o {output.tsv} &>> {log} && 
        rm {wildcards.lineage}.fixed_locus.agat.log {wildcards.lineage}.fixed_description.agat.log {wildcards.lineage}.agat.log &>> {log} || true
        """

rule fix_gff_tsv:
    input:
        tsv = rules.agat_fix_gff.output.tsv
    output:
        gff = REFDIR / "{lineage}" / "{lineage}.gff",
        tsv = REFDIR / "{lineage}" / "{lineage}.gff.tsv"
    log:
        "logs/references/fix_gff_tsv_{lineage}.log"
    shell:
        """
        python workflow/scripts/fix_gff.py -i {input.tsv} -og {output.gff} -ot {output.tsv} &> {log}
        """

