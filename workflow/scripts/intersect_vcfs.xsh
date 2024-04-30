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

    if len(vcf_files) == 1:
        print("Only one sample in lineage, no need to intersect. Copying VCF file.")
        vcf = pd.read_csv(vcf_files[0], sep='\t', comment = '#', header=None, names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample'], dtype=str)
        vcf.drop(columns=['ID', 'QUAL', 'FILTER', 'FORMAT', 'sample'], inplace=True)
        vcf.insert(vcf.columns.get_loc('POS') + 1, 'ID', '.')
        vcf.insert(vcf.columns.get_loc('ALT') + 1, 'QUAL', '.')
        vcf.insert(vcf.columns.get_loc('QUAL') + 1, 'FILTER', '.')
        vcf['var_id'] = 'var_' + lineage + '_' + (vcf.index + 1).astype(str)
        vcf['INFO'] = 'var_id=' + vcf['var_id'].astype(str) + ';' + 'MAT=1' 

        print("Creating  and saving presence table")
        sample = $(bcftools query -l @(vcf_files[0])).rstrip()
        df_presence = pd.DataFrame({'var_id': vcf['var_id'], 'sample': sample, 'lineage': lineage})
        df_presence.to_csv(presence_output, sep='\t', index=False)

        print("Saving new vcf file")
        vcf.drop(columns=['var_id'], inplace=True)
        vcf.to_csv(vcf_output, sep='\t', index=False)

    else:
        if not os.path.exists(tempdir):
            os.makedirs(tempdir)

        sites_txt_file = os.path.join(tempdir, "sites.txt")

        # Intersection
        print("Running bcftools isec")
        $(bcftools isec -p @(tempdir) @(vcf_files))

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
        df_presence['lineage'] = lineage
        df_presence.to_csv(presence_output, sep='\t', index=False)

        print("Removing temporary directory")
        $(rm -r @(tempdir))

if __name__ == '__main__':
    intersect()