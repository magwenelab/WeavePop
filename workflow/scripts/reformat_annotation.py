
import pandas as pd
import logging


log_file=snakemake.log[0]
input_tsv=snakemake.input.tsv
output_tsv=snakemake.output.tsv
output_gff=snakemake.output.gff

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading GFF table...")
    df = pd.read_csv(input_tsv, sep='\t', header=[0], dtype=str)
    
    logging.info("Defining feature number and suffix...")
    intron_mask = df['primary_tag'] == 'intron'
    df.loc[intron_mask, 'feature_number'] = df[intron_mask].groupby(['primary_tag', 'Parent']).cumcount() + 1
    df.loc[intron_mask, 'suffix'] = df.loc[intron_mask, 'primary_tag'].str[:2]
    
    logging.info("Defining new IDs...")
    df.loc[intron_mask, 'new_ID'] = df.loc[intron_mask, 'Parent'] + '-' + df.loc[intron_mask, 'suffix'] + df.loc[intron_mask, 'feature_number'].astype(int).astype(str)
    df.loc[intron_mask, 'old_ID'] = df.loc[intron_mask, 'ID']
    df.loc[intron_mask, 'ID'] = df.loc[intron_mask, 'new_ID']
    
    logging.info("Removing unnecessary columns...")
    df.drop(['feature_number', 'suffix', 'new_ID'], axis=1, inplace=True)
    
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
        "inframe_stop_codon",
        "matches_ref_protein",
        "missing_start_codon",
        "missing_stop_codon",
        "low_identity",
        "partial_mapping",
        "valid_ORF",
        "valid_ORFs"
    ]
    
    logging.info("Converting to GFF format...")
    df_gff = df.copy()
    existing_attributes = [column for column in attribute_columns if column in df.columns]

    for column in existing_attributes:
        df_gff[column] = df_gff[column].apply(lambda x: column + '=' + x if pd.notnull(x) else x)
    df_gff['attributes'] = df_gff[existing_attributes].apply(lambda x: ';'.join(x.dropna()), axis=1)
    df_gff['strand'] = df_gff['strand'].replace('-1','-').replace('1','+').replace('0','.')
    keep_columns = ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df_gff = df_gff[keep_columns]
    
    logging.info("Saving GFF...")
    df_gff.to_csv(output_gff, sep='\t', index=False, header=False)
    
    logging.info("Reformating TSV version...")
    
    logging.info("Replacing True/False for Yes/No...")
    df = df.replace({True: 'Yes', False: 'No'})
    
    logging.info("Creating start_stop_mutations column...")
    ref_mutations = ['missing_start_codon', 'missing_stop_codon', 'inframe_stop_codon']
    if any(column in df.columns for column in ref_mutations):
        exising_ref_mutations = [column for column in ref_mutations if column in df.columns]
        for column in exising_ref_mutations:
            df[column] = df[column].apply(lambda x: column if x == 'Yes' else x)
        df['start_stop_mutations'] = df[exising_ref_mutations].apply(lambda x: ','.join(x.dropna()), axis=1)
        df = df.drop(columns=exising_ref_mutations)
    
    logging.info("Renaming columns...")
    df = df.rename(columns={
        'seq_id': 'accession',
        'ID': 'feature_id',
        'Name': 'gene_name',
        'locus': 'gene_id',
        'old_ID': 'old_feature_id'})
    if 'matches_ref_protein' in df.columns:
        df.rename(columns={'matches_ref_protein': 'identical_to_main_ref'}, inplace=True)
    df.columns = df.columns.str.lower()
    
    logging.info("Reordering columns...")
    existing_columns = df.columns.tolist()
    priority_columns = [
        'accession', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame',
        'feature_id', 'gene_id', 'parent', 'gene_name',  'description', 'old_feature_id',
        'identical_to_main_ref', 'start_stop_mutations']
    
    existing_priority_columns = [column for column in priority_columns if column in existing_columns]
    
    other_columns = [column for column in existing_columns if column not in existing_priority_columns]
    df = df[existing_priority_columns + other_columns]
    
    logging.info("Final column names:")
    logging.info(df.columns.tolist())
    
    logging.info("Saving TSV...")
    df.to_csv(output_tsv, sep='\t', index=False, header=True, na_rep='')
except Exception as e:
    logging.error(f"Error: {e}")
    raise e
