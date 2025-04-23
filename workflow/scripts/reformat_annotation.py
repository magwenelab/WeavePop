import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

input_tsv=snakemake.input.tsv
output_tsv=snakemake.output.tsv
output_gff=snakemake.output.gff
lineage=snakemake.params.lineage

print("Reading GFF table...")
df = pd.read_csv(input_tsv, sep='\t', header=[0], dtype=str)
print("Sorting the accessions because AGAT gff2tsv disorders them...")
df.sort_values(by=["seq_id"], kind='mergesort', inplace=True)

print("Defining feature number and suffix...")
intron_mask = df['primary_tag'] == 'intron'
df.loc[intron_mask, 'feature_number'] = df[intron_mask].groupby(['primary_tag', 'Parent']).cumcount() + 1
df.loc[intron_mask, 'suffix'] = df.loc[intron_mask, 'primary_tag'].str[:2]

print("Defining new IDs...")
df.loc[intron_mask, 'new_ID'] = df.loc[intron_mask, 'Parent'] + '-' + df.loc[intron_mask, 'suffix'] + df.loc[intron_mask, 'feature_number'].astype(int).astype(str)
df.loc[intron_mask, 'old_ID'] = df.loc[intron_mask, 'ID']
df.loc[intron_mask, 'ID'] = df.loc[intron_mask, 'new_ID']

print("Removing unnecessary columns...")
df.drop(['feature_number', 'suffix', 'new_ID'], axis=1, inplace=True)

print("Replacing True/False for Yes/No...")
df = df.replace({'True': 'Yes', 'False': 'No'})

print("Creating start_stop_mutations column...")
ref_mutations = ['missing_start_codon', 'missing_stop_codon', 'inframe_stop_codon']
if any(column in df.columns for column in ref_mutations):
    exising_ref_mutations = [column for column in ref_mutations if column in df.columns]
    for column in exising_ref_mutations:
        df[column] = df[column].apply(lambda x: column if x == 'Yes' else x)
    df['start_stop_mutations'] = df[exising_ref_mutations].apply(lambda x: ','.join(x.dropna()), axis=1)
    df['start_stop_mutations'] = df['start_stop_mutations'].replace('', None)
    print(set(df['start_stop_mutations']))
    df = df.drop(columns=exising_ref_mutations)
    
print("Renaming matches_ref_protein to indentical_to_main_ref if present...")
if 'matches_ref_protein' in df.columns:
    df.rename(columns={'matches_ref_protein': 'identical_to_main_ref'}, inplace=True)
    
print("Converting to GFF format...")
attribute_columns = [
"ID",
"locus",
"Parent",
"Name",
"description",
"old_ID",
"sequence_ID",
"copy_num_ID",
"coverage",
"extra_copy_number",
"identical_to_main_ref",
"start_stop_mutations",
"low_identity",
"partial_mapping",
"valid_ORF",
"valid_ORFs",
"repeat_fraction"]

existing_attributes = [column for column in attribute_columns if column in df.columns]

df_gff = df.copy()

for column in existing_attributes:
    df_gff[column] = df_gff[column].apply(lambda x: column + '=' + x if pd.notnull(x) else x)
df_gff['attributes'] = df_gff[existing_attributes].apply(lambda x: ';'.join(x.dropna()), axis=1)
df_gff['strand'] = df_gff['strand'].replace('-1','-').replace('1','+').replace('0','.')
keep_columns = ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
df_gff = df_gff[keep_columns]

print("Saving GFF...")
df_gff.to_csv(output_gff, sep='\t', index=False, header=False)

print("Renaming columns of TSV version...")
df = df.rename(columns={
    'seq_id': 'accession',
    'ID': 'feature_id',
    'Name': 'gene_name',
    'locus': 'gene_id',
    'old_ID': 'old_feature_id'})

df.columns = df.columns.str.lower()

print("Reordering columns...")
existing_columns = df.columns.tolist()
priority_columns = [
    'accession', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame',
    'feature_id', 'gene_id', 'parent', 'gene_name',  'description', 'old_feature_id',
    'identical_to_main_ref', 'start_stop_mutations']

existing_priority_columns = [column for column in priority_columns if column in existing_columns]

other_columns = [column for column in existing_columns if column not in existing_priority_columns]
df = df[existing_priority_columns + other_columns]

print("Add new column with lineage name...")
df['lineage'] = lineage

print("Final column names:")
print(df.columns.tolist())

print("Saving TSV...")
df.to_csv(output_tsv, sep='\t', index=False, header=True, na_rep='')
