import click
import pandas as pd
import io
import os

@click.command()
@click.option('--gff_file', '-g', help='Input GFF file of the reference genome', required=True, type=click.Path(exists=True))
@click.option('--lineage', '-l', help='Lineage', required=True, type=str)
@click.option('--output_gff', '-og', help='Path for GFF output file with variants', required=True, type=click.Path())
@click.option('--output_vars', '-ov', help='Path for variant output file with genes', required=True, type=click.Path())
@click.option('--temp_dir', '-t', help='Path of the temporary directory', default='temp/', type=click.Path())
@click.argument('vcf_files', nargs=-1, type=click.Path(exists=True), required=True)


def variant_intersection(gff_file, lineage, vcf_files, output_gff, output_vars, temp_dir):
    print("Obtaining intersection of variants from the VCF files")
    temp = os.path.join(temp_dir, "bcf_isec_" + lineage)
    
    $(bcftools isec -p @(temp) @(vcf_files))

    file_path = temp + '/sites.txt'
    var_original = pd.read_csv(file_path, sep='\t', header=None, names=['CHROM', 'POS', 'REF', 'ALT', 'PRESENCE'], dtype=str)
    
    sample_names = []
    for vcf_file in vcf_files:
        sample_names.append($(bcftools query -l @(vcf_file)).rstrip())

    presence_columns = var_original['PRESENCE'].apply(lambda x: pd.Series(list(x)))
    presence_columns.columns = sample_names
    var_presence = pd.concat([var_original, presence_columns], axis=1)
    var_presence = var_presence.drop(columns=['PRESENCE'])
    var_presence['VAR'] = 'Var' + var_presence.index.astype(str)

    var_bed = var_presence[['CHROM', 'POS', 'POS', 'VAR']].copy()
    var_bed.columns = ['CHROM', 'POS', 'POS2', 'VAR']
    temp_bed_path = os.path.join(temp_dir, 'var_bed.tsv')
    var_bed.to_csv(temp_bed_path, sep='\t', index=False, header=False)

    print("Preparing GFF file for intersection")
    gff_df = pd.read_csv(gff_file, sep='\t', header=0, comment='#', dtype=str)  
    if 'ID' not in gff_df.columns or not gff_df['ID'].is_unique:
        gff_df['unique_id'] = range(1, len(gff_df) + 1)
        gff_df['unique_id'] = 'ID_' + gff_df['unique_id'].astype(str)
    else:
        gff_df['unique_id'] = gff_df['ID']

    column_names = ['seq_id', 'start', 'end', 'unique_id']
    existing_columns = [col for col in column_names if col in gff_df.columns]
    gff_temp = gff_df[existing_columns]
    temp_gff_path = os.path.join(temp_dir, 'gff_df.tsv')
    gff_temp.to_csv(temp_gff_path, sep='\t', index=False, header=False)

    print("Intersecting variants with GFF file")
    intersection = $(bedtools intersect -wa -wb -a @(temp_bed_path) -b @(temp_gff_path))
    intersection_df = pd.read_csv(io.StringIO(intersection), sep='\t', header=None, dtype=str)
    intersect_columns = var_bed.columns.tolist() + gff_temp.columns.tolist()
    intersection_df.columns = intersect_columns
    intersection_df = intersection_df.drop(columns=['CHROM', 'POS','POS2', 'seq_id', 'start', 'end'])

    print("Merging GFF file with variants")
    feature_var = intersection_df.groupby('unique_id').agg(lambda x: ';'.join(x.unique()))
    gff_var = gff_df.merge(feature_var, on='unique_id', how='left')
    
    print("Saving GFF file with variants")
    gff_var.to_csv(output_gff, sep='\t', index=False, quoting=1, quotechar='"')

    print("Merging variant presence with gene IDs")
    column_names = ['primary_tag','unique_id', 'Name', 'description', 'product', 'locus_tag']
    existing_columns = [col for col in column_names if col in gff_df.columns]
    gene_ID = gff_df[existing_columns]
    gene_ID = gene_ID[gene_ID['primary_tag'].isin(['protein_coding_gene', 'gene', 'pseudogene', 'ncRNA_gene'])]
    var_gene = gene_ID.merge(intersection_df, on='unique_id', how='left')
    variant_presence_gene = var_presence.merge(var_gene, on='VAR', how='left')

    print("Saving variant presence with gene IDs")
    variant_presence_gene.to_csv(output_vars, sep='\t', index=False, quoting=1, quotechar='"')

    print("Removing temporary files")
    $(rm -r @(temp_dir))

if __name__ == '__main__':
    variant_intersection()