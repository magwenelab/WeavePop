from pathlib import Path
import pandas as pd
import logging

log_file=snakemake.log[0]
input=snakemake.input
output=snakemake.output[0]

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Create dictionary with lineage and path...")
    paths = {}
    for path in input:
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
        'accession', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 
        'feature_id', 'parent', 'gene_id', 'gene_name', 'description', 'old_feature_id', 'repeat_fraction', 'lineage',
        'identical_to_main_ref','start_stop_mutations']
    existing_columns = [column for column in keep_columns if column in df.columns]
    logging.info("Keeping columns:")
    logging.info(existing_columns)
    df = df[existing_columns]

    logging.info("Saving...")
    df.to_csv(output, sep='\t', index=False)
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e
