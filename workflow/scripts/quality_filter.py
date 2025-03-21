
import pandas as pd
import logging
import os
import sys

log_file=snakemake.log[0]

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

input_stats=snakemake.input[0]
input_chromosomes=snakemake.input[1]
filter=snakemake.params.filter
output_stats=snakemake.output.stats
output_metadata=snakemake.output.metadata
output_chromosomes=snakemake.output.chromosomes
metadata = snakemake.params.metadata


try:
    logging.info("Reading tables...")
    stats = pd.read_csv(input_stats, sep="\t", header=0)
    chromosomes = pd.read_csv(input_chromosomes, header=0)
    if filter:
        logging.info("Filtering samples...")
        stats_filtered = stats[stats["quality_warning"].isna()]
        metadata_filtered = metadata.loc[metadata["sample"].isin(stats_filtered["sample"]),]
        chromosomes_filtered = chromosomes.loc[chromosomes["lineage"].isin(metadata_filtered["lineage"]),]
        logging.info("Successfully filtered samples from tables.\n")
    else:
        logging.info("No filtering requested, copying tables...")
        stats_filtered = stats
        metadata_filtered = metadata
        chromosomes_filtered = chromosomes
    
    logging.info("Cleaning column names...")
    metadata_filtered.columns = metadata_filtered.columns.str.lower()
    metadata_filtered.columns = metadata_filtered.columns.str.replace(' ', '_')
    logging.info("Adding dataset column...")
    metadata_filtered.loc[:, 'dataset'] = "X"
        
    logging.info("Saving filtered tables...")
    stats_filtered.to_csv(output_stats, index=False, header=True, sep = "\t")
    metadata_filtered.to_csv(output_metadata, index=False)
    chromosomes_filtered.to_csv(output_chromosomes, index=False)
    logging.info("Done!")
except Exception as e:
    log_file.write(f"Error: {e}\n")
    raise e