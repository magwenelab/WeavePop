import pandas as pd
import vcf
import click
import io
import os
import duckdb

# Variables
print("Setting up variables")
lineage = 'VNBI'
temp_dir = 'temp_VNBI'
vcf_files=('results/samples/snippy/SRS8318900/snps.vcf.gz','results/samples/snippy/SRS8318901/snps.vcf.gz')
db_name = 'Cryptococcus_neoformans_VNBI'

temp = os.path.join(temp_dir, "bcf_isec_" + lineage)
# Make temp directory if it is not present
if not os.path.exists(temp):
    os.makedirs(temp)

sites_txt_file = temp + '/sites.txt'
sites_vcf_path = temp + '/sites.vcf'
ann_vcf_path = os.path.join(temp_dir,'snps.ann.vcf')
bcftools_log_path = os.path.join(temp_dir, 'bcftools.log')
snpeff_log_path = os.path.join(temp_dir, 'snpeff.log') 
snpeff_html_path = os.path.join(temp_dir, 'snpeff.html')
presence_path = os.path.join(temp_dir, 'presence.csv')
variants_path = os.path.join(temp_dir, 'variants.csv')
effects_path = os.path.join(temp_dir, 'effects.csv')
lofs_path = os.path.join(temp_dir, 'lofs.csv')
nmds_path = os.path.join(temp_dir, 'nmds.csv')
genes_path = os.path.join(temp_dir, 'genes.csv')

# Intersection
print("Running bcftools isec")
$(bcftools isec -p @(temp) @(vcf_files) 2> @(bcftools_log_path))

# Convert sites.txt to VCF
print("Converting sites.txt to VCF")
sites_txt = pd.read_csv(sites_txt_file, sep='\t', header=None, names=['CHROM', 'POS', 'REF', 'ALT', 'INFO'], dtype=str)
sites_txt['VAR_ID'] = 'VAR_' + (sites_txt.index + 1).astype(str)
sites_vcf = sites_txt.copy()
sites_vcf['INFO'] = 'VAR_ID=' + sites_vcf['VAR_ID'].astype(str) + ';' + 'MAT=' + sites_vcf['INFO'].astype(str)
sites_vcf.drop(columns=['VAR_ID'], inplace=True)
sites_vcf.insert(sites_vcf.columns.get_loc('POS') + 1, 'ID', '.')
sites_vcf.insert(sites_vcf.columns.get_loc('ALT') + 1, 'QUAL', '.')
sites_vcf.insert(sites_vcf.columns.get_loc('QUAL') + 1, 'FILTER', '.')
sites_vcf.to_csv(sites_vcf_path, sep='\t', index=False)

#Get presence/absence of variants
print("Getting presence/absence of variants data frame")
sample_names = []
for vcf_file in vcf_files:
    sample_names.append($(bcftools query -l @(vcf_file) 2>> @(bcftools_log_path)).rstrip())

presence_columns = sites_txt['INFO'].apply(lambda x: pd.Series(list(x)))
presence_columns.columns = sample_names
presence_matrix = pd.concat([sites_txt, presence_columns], axis=1)
presence_matrix = presence_matrix.drop(columns=['INFO', 'CHROM', 'POS', 'REF', 'ALT'])
# Transpose df_presence
presence_melt = presence_matrix.melt(id_vars='VAR_ID', var_name='sample', value_name='value')
df_presence = presence_melt[presence_melt['value'] == '1'].copy()
df_presence.drop(columns='value', inplace=True)

# Run SnpEFF
print("Running SnpEff")
$(snpEff ann -v -classic -s @(snpeff_html_path) @(db_name) @(sites_vcf_path) 1> @(ann_vcf_path) 2> @(snpeff_log_path))

# Get tables from SnpEFF result
print("Getting tables from SnpEff result")
# Create a list to store the variant data
data_effects = []
data_variants = []
data_lofs = []
data_nmds = []
# Iterate over the records in the VCF file
with open(ann_vcf_path, 'r') as file:
    # Open the VCF file
    vcf_reader = vcf.Reader(file)

    # Add a unique identifier to each record
    for i, record in enumerate(vcf_reader):
        var_id = record.INFO['VAR_ID']
        var_id = var_id[0]
        data_variants.append([var_id, record.CHROM, record.POS, record.REF, record.ALT])
                
        # Get the SnpEff annotations from the INFO field
        annotations = record.INFO['EFF']
        # Get the LOF annotations from the INFO field
        lofs = record.INFO.get('LOF', [])
        # Get the NMD annotations from the INFO field
        nmds = record.INFO.get('NMD', [])
        # Iterate over all annotations
        for annotation in annotations:
            # Split the annotation into its components
            annotation_parts = annotation.split('(')
            
            # Get the gene name and variant effect
            effect_type = annotation_parts[0]
            effect_info = annotation_parts[1]
            effect_info_parts = effect_info.split('|')
            impact = effect_info_parts[0]
            effect = effect_info_parts[1]
            codon_change = effect_info_parts[2]
            aa_change = effect_info_parts[3]
            aa_length = effect_info_parts[4]
            gene_name = effect_info_parts[5]
            biotype = effect_info_parts[6]
            coding = effect_info_parts[7]
            transcript_id = effect_info_parts[8]
            exon_rank = effect_info_parts[9]
            data_effects.append([var_id,  effect_type, impact, effect, codon_change, aa_change, aa_length, gene_name, biotype, coding, transcript_id, exon_rank])
        # Iterate over all lofs
        for lof in lofs:
            lof_parts = lof.split('|')
            first = lof_parts[0]
            gene_name = first.replace('(', '')
            gene_ID = lof_parts[1]
            nb_transcripts = lof_parts[2]
            last = lof_parts[3]
            percent_transcripts = last.replace(')', '')
            data_lofs.append([var_id, gene_name, gene_ID, nb_transcripts, percent_transcripts])
        # Iterate over all nmds
        for nmd in nmds:
            nmd_parts = nmd.split('|')
            first = nmd_parts[0]
            gene_name = first.replace('(', '')
            gene_ID = nmd_parts[1]
            nb_transcripts = nmd_parts[2]
            last = nmd_parts[3]
            percent_transcripts = last.replace(')', '')
            data_nmds.append([var_id, gene_name, gene_ID, nb_transcripts, percent_transcripts])

# Create dataframes from the data objects
print("Creating dataframes")
df_variants = pd.DataFrame(data_variants, columns=['VAR_ID', 'Chrom', 'Pos', 'Ref', 'Alt'])

df_effects = pd.DataFrame(data_effects, columns=['VAR_ID', 'Type', 'Impact','Effect', 'Codon_change', 'Amino_acid_change', 'Amino_acid_length', 'Gene_name', 'Transcript_biotype', 'Gene_coding', 'Transcript_ID', 'Exon_rank'])
df_effects['Effect_ID'] = 'EFF_' + (df_effects.index + 1).astype(str)
df_effects.insert(0, 'Effect_ID', df_effects.pop('Effect_ID'))

df_lofs = pd.DataFrame(data_lofs, columns=['VAR_ID', 'Gene_name','Gene_ID', 'Nb_transcripts', 'Percent_transcripts'])

df_nmds = pd.DataFrame(data_nmds, columns=['VAR_ID', 'Gene_name','Gene_ID', 'Nb_transcripts', 'Percent_transcripts'])

df_genes = df_effects[['Gene_name','Transcript_ID', 'Effect_ID']].copy()
df_genes.replace('', None , inplace=True)
df_genes.dropna(axis = 0, inplace=True)
df_genes.sort_values('Transcript_ID', inplace=True)

# Save the tables to CSV files
print("Saving dataframes to CSV files")
df_presence.to_csv(presence_path, index=False)
df_variants.to_csv(variants_path, index=False)
df_effects.to_csv(effects_path, index=False)
df_lofs.to_csv(lofs_path, index=False)
df_nmds.to_csv(nmds_path, index=False)
df_genes.to_csv(genes_path, index=False)

# DUCKDB

# Connect to the database
con = duckdb.connect(database='temp_VNBI/duck_db_VNBI.db')
con.execute("CREATE TABLE IF NOT EXISTS presence AS SELECT * FROM df_presence")
con.execute("CREATE TABLE IF NOT EXISTS variants AS SELECT * FROM df_variants")
con.execute("CREATE TABLE IF NOT EXISTS effects AS SELECT * FROM df_effects")
con.execute("CREATE TABLE IF NOT EXISTS lofs AS SELECT * FROM df_lofs")
con.execute("CREATE TABLE IF NOT EXISTS nmds AS SELECT * FROM df_nmds")
con.execute("CREATE TABLE IF NOT EXISTS genes AS SELECT * FROM df_genes")
