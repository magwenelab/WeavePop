import pandas as pd
import click
import io
import os
import re

@click.command()
@click.option('--vcf_output', '-v', type = click.Path(), help='Path to output intersected VCF file')
@click.option('--presence_output', '-p', type = click.Path(), help='Path to output variant presence table')
@click.option('--lineage', '-l', type = str, help='Lineage to intersect')
@click.option('--tempdir', '-t', type = click.Path(), help='Output directory for intermediate files')
@click.argument('vcf_files', nargs = -1, type=click.Path())
def intersect(vcf_output, presence_output, lineage, tempdir, vcf_files):
    dir = os.path.join(tempdir, str("bcf_isec_" + lineage))
    if not os.path.exists(dir):
        os.makedirs(dir)

    sites_txt_file = os.path.join(dir, "sites.txt")

    # Intersection
    print("Running bcftools isec")
    $(bcftools isec -p @(dir) @(vcf_files))

    # Convert sites.txt to VCF
    print("Converting sites.txt to VCF")
    sites_txt = pd.read_csv(sites_txt_file, sep='\t', header=None, names=['#CHROM', 'POS', 'REF', 'ALT', 'INFO'], dtype=str)
    sites_txt['var_id'] = 'var_' + lineage + '_' + (sites_txt.index + 1).astype(str)
    sites_vcf = sites_txt.copy()
    sites_vcf['INFO'] = 'var_id=' + sites_vcf['var_id'].astype(str) + ';' + 'MAT=' + sites_vcf['INFO'].astype(str)
    sites_vcf.drop(columns=['var_id'], inplace=True)
    sites_vcf.insert(sites_vcf.columns.get_loc('POS') + 1, 'ID', '.')
    sites_vcf.insert(sites_vcf.columns.get_loc('ALT') + 1, 'QUAL', '.')
    sites_vcf.insert(sites_vcf.columns.get_loc('QUAL') + 1, 'FILTER', '.')
    sites_vcf.to_csv(vcf_output, sep='\t', index=False)

    #Get presence/absence of variants
    print("Adding sample names to presence column")
    sample_names = []
    for vcf_file in vcf_files:
        sample_names.append($(bcftools query -l @(vcf_file)).rstrip())
    print(sample_names)
    print("Creating presence columns")
    presence_columns = sites_txt['INFO'].apply(lambda x: pd.Series(list(x)))
    presence_columns.columns = sample_names
    print("Creating presence matrix")
    presence_matrix = pd.concat([sites_txt, presence_columns], axis=1)
    presence_matrix = presence_matrix.drop(columns=['INFO', '#CHROM', 'POS', 'REF', 'ALT'])
    print("Converting matrix to dataframe")
    presence_melt = presence_matrix.melt(id_vars='var_id', var_name='sample', value_name='value')
    df_presence = presence_melt[presence_melt['value'] == '1'].copy()
    df_presence.drop(columns='value', inplace=True)
    df_presence.to_csv(presence_output, sep='\t', index=False)

if __name__ == '__main__':
    intersect()