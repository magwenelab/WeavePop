import click
import pandas as pd
import io
import os

@click.command()
@click.option('--gff-file', '-g', help='Input GFF file of the reference genome')
@click.option('--lineage', '-l', help='Lineage', required=True, type=str)
@click.argument('vcf_files', nargs=-1, type=click.Path(exists=True), required=True)

def variant_intersection(gff_file, lineage, vcf_files):
    # Obtain a presence/absence matrix of variants for each sample and a file with the variant information
    $(bcftools isec -p bcf_isec_@(lineage) @(vcf_files))

    file_path = os.path.join('bcf_isec_' + lineage, 'sites.txt')
    df = pd.read_csv(file_path, sep='\t', header=None, names=['CHROM', 'POS', 'REF', 'ALT', 'PRESENCE'], dtype=str)
    
    sample_names = []
    for vcf_file in vcf_files:
        sample_names.append($(bcftools query -l @(vcf_file)).rstrip())

    presence_columns = df['PRESENCE'].apply(lambda x: pd.Series(list(x)))
    presence_columns.columns = sample_names
    df = pd.concat([df, presence_columns], axis=1)
    df = df.drop(columns=['PRESENCE'])
    df['VAR'] = 'Var' + df.index.astype(str)

    var_bed = df[['CHROM', 'POS', 'POS']].copy()
    var_bed.columns = ['CHROM', 'POS', 'POS2']
    var_bed.to_csv('var_bed.tsv', sep='\t', index=False, header=False)

    # Obtain the gene information for each variant
    gff_df = pd.read_csv(gff_file, sep='\t', header=0, comment='#', dtype=str)
    gff_df = gff_df[(gff_df['primary_tag'] == 'gene') | (gff_df['primary_tag'] == 'protein_coding_gene')]
    
    column_names = ['seq_id', 'start', 'end', 'ID', 'Name', 'description', 'product', 'locus_tag']
    existing_columns = [col for col in column_names if col in gff_df.columns]
    gff_df = gff_df[existing_columns]
    gff_df.to_csv('gff_df.tsv', sep='\t', index=False, header=False)

    intersect_columns = ['VAR'] + gff_df.columns.tolist() + var_bed.columns.tolist() + sample_names

    gene_variant = $(bedtools intersect -wa -wb -a var_bed.tsv -b gff_df.tsv)     
    gene_variant_df = pd.read_csv(io.StringIO(gene_variant), sep='\t', header=None, dtype=str)
    gene_variant_df.columns = intersect_columns
    gene_variant_df = gene_variant_df.drop(columns=['POS2', 'seq_id', 'start', 'end'])

    merged_df = df.merge(gene_variant_df, on=['CHROM', 'POS'], how='left')
    merged_df = merged_df[['VAR'] + [col for col in merged_df.columns if col != 'VAR']]

    sample_columns = [col for col in merged_df.columns if col in sample_names]
    other_columns = [col for col in merged_df.columns if col not in sample_names]

    merged_df = merged_df[other_columns + sample_columns]

    merged_df.to_csv('merged_df.tsv', sep='\t', index=False)

    if __name__ == '__main__':
    variant_intersection()

    # Its missing to set the output file names correctly and the tenmporary files
    