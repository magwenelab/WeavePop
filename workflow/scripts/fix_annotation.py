
import pandas as pd
from pathlib import Path
import logging

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

log_file=snakemake.log[0]
input_tsv=snakemake.input.tsv
output_tsv=snakemake.output.tsv
output_gff=snakemake.output.gff

try:
    logging.info("Reading GFF table...")
    df=pd.read_csv(input_tsv, sep='\t', header=[0])
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
    logging.info("Saving TSV...")
    df.to_csv(output_tsv, sep='\t', index=False, header=True, na_rep='')
    
    logging.info("Making GFF...")
    attribute_columns = [
        "copy_num_ID",
        "coverage",
        "description",
        "extra_copy_number",
        "ID",
        "inframe_stop_codon",
        "locus",
        "matches_ref_protein",
        "missing_start_codon",
        "missing_stop_codon",
        "Name",
        "old_ID",
        "Parent",
        "sequence_ID",
        "valid_ORF",
        "valid_ORFs"
    ]
    existing_attributes = [column for column in attribute_columns if column in df.columns]
    df_gff = df.copy()
    for column in existing_attributes:
        df_gff[column] = df_gff[column].apply(lambda x: column + '=' + x if pd.notnull(x) else x)
    df_gff['attributes'] = df_gff[existing_attributes].apply(lambda x: ';'.join(x.dropna()), axis=1)
    df_gff.drop(existing_attributes, axis=1, inplace=True)
    df_gff['strand'] = df_gff['strand'].replace(-1,'-').replace(1,'+')
    logging.info("Saving GFF...")
    df_gff.to_csv(output_gff, sep='\t', index=False, header=False)
except Exception as e:
    logging.error(f"Error: {e}")
    raise e
