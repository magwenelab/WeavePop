import pandas as pd
import click
import io
import os
import re

@click.command()
@click.option('--lineage', type = str, help='Lineage to intersect')
@click.option('--outdir', type = click.Path(), help='Path to output directory')
@click.argument('vcf_files', help='List of VCF files to intersect')

def intersect(lineage, vcf_files, outdir)
    temp = os.path.join(temp_dir, "bcf_isec_" + lineage)
    if not os.path.exists(temp):
        os.makedirs(temp)

    sites_txt_file = temp + '/sites.txt'
    sites_vcf_path = temp + '/sites.vcf'
    ann_vcf_path = os.path.join(temp_dir,'snps.ann.vcf')
    bcftools_log_path = os.path.join(temp_dir, 'bcftools.log')
    snpeff_log_path = os.path.join(temp_dir, 'snpeff.log') 
    snpeff_html_path = os.path.join(temp_dir, 'snpeff.html')

    # Intersection
    print("Running bcftools isec")
    $(bcftools isec -p @(temp) @(vcf_files) 2> @(bcftools_log_path))

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
    sites_vcf.to_csv(sites_vcf_path, sep='\t', index=False)

    #Get presence/absence of variants
    print("Adding sample names to presence column")
    sample_names = []
    for vcf_file in vcf_files:
        sample_names.append($(bcftools query -l @(vcf_file) 2>> @(bcftools_log_path)).rstrip())
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