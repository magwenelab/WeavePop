import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import pandas as pd
import click
import io
import os
import re

vcf_files=snakemake.input.vcfs

vcf_output=snakemake.output.vcf
presence_output=snakemake.output.tsv

lineage=snakemake.wildcards.lineage
tempdir=snakemake.params.tmp_dir

print("Creating comments for the top of the VCF file...")
software_legend="##WeavePop_Script=intersect_vcfs.xsh"
var_id_legend=("##INFO=<ID=var_id,Number=1,Type=String,Description="+
    "\"Identifier of unique lineage variants.\">")


if len(vcf_files) == 1:
    print("Only one sample in lineage, no need to intersect. Copying VCF file...")
    vcf = pd.read_csv(
        vcf_files[0], sep='\t', comment='#', header=None,
        names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER','INFO', 'FORMAT', 'sample'],
        dtype=str)
    vcf.drop(columns=['ID', 'QUAL', 'FILTER', 'FORMAT', 'sample'], inplace=True)
    vcf.insert(vcf.columns.get_loc('POS') + 1, 'ID', '.')
    vcf.insert(vcf.columns.get_loc('ALT') + 1, 'QUAL', '.')
    vcf.insert(vcf.columns.get_loc('QUAL') + 1, 'FILTER', '.')
    vcf['var_id'] = 'var_' + lineage + '_' + (vcf.index + 1).astype(str)
    vcf['INFO'] = 'var_id=' + vcf['var_id'].astype(str) + ';' + 'MAT=1' 

    print("Obtaining sample name from the VCF file...")
    sample = $(bcftools query -l @(vcf_files[0]) 2>> @(log_file)).rstrip()
    sample_names_legend=("##INFO=<ID=MAT,Number=.,Type=String,Description=" +
        "\"Presence/Absence matrix for samples. Order:" +
        sample +
        '">')

    print("Creating  and saving presence table...")
    df_presence = pd.DataFrame({'var_id': vcf['var_id'], 'sample': sample})
    df_presence.to_csv(presence_output, sep='\t', index=False)

    print("Saving new vcf file...")
    vcf.drop(columns=['var_id'], inplace=True)
    vcf.to_csv(vcf_output, sep='\t', index=False)

else:
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    sites_txt_file = os.path.join(tempdir, "sites.txt")

    print("Running bcftools isec to intersect the variants of all the samples in the lineage...")
    $(bcftools isec -p @(tempdir) @(vcf_files) 2>> @(log_file))

    print("Converting sites.txt to VCF...")
    sites_txt = pd.read_csv(sites_txt_file, sep='\t', 
        header=None,
        names=['#CHROM', 'POS', 'REF', 'ALT', 'INFO'],
        dtype=str)
    sites_txt['var_id'] = 'var_' + lineage + '_' + (sites_txt.index + 1).astype(str)
    sites_vcf = sites_txt.copy()
    sites_vcf['INFO'] = (
        'var_id=' + sites_vcf['var_id'].astype(str) +
        ';' +
        'MAT=' + sites_vcf['INFO'].astype(str))
    sites_vcf.drop(columns=['var_id'], inplace=True)
    sites_vcf.insert(sites_vcf.columns.get_loc('POS') + 1, 'ID', '.')
    sites_vcf.insert(sites_vcf.columns.get_loc('ALT') + 1, 'QUAL', '.')
    sites_vcf.insert(sites_vcf.columns.get_loc('QUAL') + 1, 'FILTER', '.')

    print("Saving the intersected VCF file...")
    sites_vcf.to_csv(vcf_output, sep='\t', index=False)

    print("Obtaining sample names from the order of the original VCF files...")
    sample_names = []
    for vcf_file in vcf_files:
        sample_names.append($(bcftools query -l @(vcf_file) 2>> @(log_file)).rstrip())
    sample_names_legend=("##INFO=<ID=MAT,Number=.,Type=String,Description=" +
        "\"Presence/Absence matrix for samples. Order:" +
        ','.join(sample_names) +
        '">')

    print("Creating presence columns from the sites.txt file...")
    presence_columns = sites_txt['INFO'].apply(lambda x: pd.Series(list(x)))
    presence_columns.columns = sample_names
    print("Creating presence matrix...")
    presence_matrix = pd.concat([sites_txt, presence_columns], axis=1)
    presence_matrix = presence_matrix.drop(columns=['INFO', '#CHROM', 'POS', 'REF', 'ALT'])
    print("Converting matrix to dataframe...")
    presence_melt = presence_matrix.melt(id_vars='var_id', var_name='sample', value_name='value')
    df_presence = presence_melt[presence_melt['value'] == '1'].copy()
    df_presence.drop(columns='value', inplace=True)
    print("Saving presence table...")
    df_presence.to_csv(presence_output, sep='\t', index=False)

    print("Removing temporary directory...")
    $(rm -r @(tempdir))

print("Adding comments to the top of the VCF file...")
with open(vcf_output, 'r') as vcf_read_file:
    vcf_content = vcf_read_file.readlines()

vcf_content.insert(0, software_legend + '\n')
vcf_content.insert(1, var_id_legend + '\n')
vcf_content.insert(2, sample_names_legend + '\n')

with open(vcf_output, 'w') as vcf_write_file:
    vcf_write_file.writelines(vcf_content)  

print("Done!")