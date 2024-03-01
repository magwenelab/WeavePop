import click
import pandas as pd
import io
import os

@click.command()
@click.option('--gff_file', '-g', help='Input GFF file of the reference genome', required=True, type=click.Path(exists=True))
@click.option('--lineage', '-l', help='Lineage', required=True, type=str)
@click.option('--output', '-o', help='Output file path', required=True, type=click.Path())
@click.option('--temp_dir', '-t', help='Path of the temporary directory', default='temp/', type=click.Path())
@click.argument('vcf_files', nargs=-1, type=click.Path(exists=True), required=True)


def variant_intersection(gff_file, lineage, vcf_files, output, temp_dir):
    print("Obtaining intersection of variants from the VCF files")
    temp = os.path.join(temp_dir, "bcf_isec_" + lineage)
    
    $(bcftools isec -p @(temp) @(vcf_files))

    file_path = temp + '/sites.txt'
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
    temp_bed = os.path.join(temp_dir, 'var_bed.tsv')
    var_bed.to_csv(temp_bed, sep='\t', index=False, header=False)

    print("Obtaining genes where variants are located")
    gff_df = pd.read_csv(gff_file, sep='\t', header=0, comment='#', dtype=str)
    gff_df = gff_df[(gff_df['primary_tag'] == 'gene') | (gff_df['primary_tag'] == 'protein_coding_gene')]
    
    column_names = ['seq_id', 'start', 'end', 'ID', 'Name', 'description', 'product', 'locus_tag']
    existing_columns = [col for col in column_names if col in gff_df.columns]
    gff_df = gff_df[existing_columns]
    temp_gff = os.path.join(temp_dir, 'gff_df.tsv')
    gff_df.to_csv(temp_gff, sep='\t', index=False, header=False)

    print("Intersecting variants with genes")
    gene_variant = $(bedtools intersect -wa -wb -a @(temp_bed) -b @(temp_gff))
    gene_variant_df = pd.read_csv(io.StringIO(gene_variant), sep='\t', header=None, dtype=str)
    intersect_columns = var_bed.columns.tolist() + gff_df.columns.tolist()
    gene_variant_df.columns = intersect_columns
    gene_variant_df = gene_variant_df.drop(columns=['POS2', 'seq_id', 'start', 'end'])

    print("Creating output file")
    merged_df = df.merge(gene_variant_df, on=['CHROM', 'POS'], how='left')
    merged_df = merged_df[['VAR'] + [col for col in merged_df.columns if col != 'VAR']]
    sample_columns = [col for col in merged_df.columns if col in sample_names]
    other_columns = [col for col in merged_df.columns if col not in sample_names]
    merged_df = merged_df[other_columns + sample_columns]

    print("Saving output file")
    merged_df.to_csv(output, sep='\t', index=False)

    print("Removing temporary files")
    $(rm -r @(temp_dir))

if __name__ == '__main__':
    variant_intersection()

