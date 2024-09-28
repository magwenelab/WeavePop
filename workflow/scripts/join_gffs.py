from pathlib import Path
import pandas as pd
import logging

log_file = snakemake.log[0]

# Configure logging
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Create dictionary with lineage and path...")
    paths = {}
    for path in snakemake.input:
        lineage = Path(Path(path).stem).stem
        paths[lineage] = path   
        logging.info("Available lineages and paths:") 
        logging.info(paths)
    logging.info("Files to concatenate:")    
    logging.info(paths)
    logging.info("Read files and concatenate dataframes into list...")
    dfs = []
    for lineage, path in paths.items():
        df = pd.read_csv(path, sep='\t', low_memory=False, header=0)
        df['lineage'] = lineage
        dfs.append(df)
    logging.info("Concatenate dataframes into one...")
    df = pd.concat(dfs)
    keep_columns = [
        'seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 
        'frame', 'ID', 'Parent', 'locus', 'Name', 'description', 'old_ID', 'lineage',
        'matches_ref_protein','missing_start_codon', 'missing_stop_codon', 'inframe_stop_codon']
    existing_columns = [column for column in keep_columns if column in df.columns]
    logging.info("Keeping columns:")
    logging.info(existing_columns)
    df = df[existing_columns]
    logging.info("Replacing True/False for Yes/No...")
    df = df.replace({True: 'Yes', False: 'No'})
    logging.info("Renaming columns...")
    df = df.rename(columns={
        'seq_id': 'accession',
        'ID': 'feature_id',
        'Name': 'gene_name',
        'locus': 'gene_id',
        'old_ID': 'old_feature_id'})
    logging.info("Saving...")
    df.to_csv(snakemake.output[0], sep='\t', index=False)
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e
