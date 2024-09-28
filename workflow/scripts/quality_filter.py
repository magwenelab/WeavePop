
import pandas as pd
import logging

log_file = snakemake.log[0]

# Configure logging
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading tables...")
    stats = pd.read_csv(snakemake.input[0], sep="\t", header = 0)
    metadata = pd.read_csv(snakemake.input[1], header=0)
    if snakemake.params.exclude:
        logging.info("Filtering samples...")
        stats_filtered = stats[stats["quality_warning"].isna()]
        metadata_filtered = metadata.loc[metadata["sample"].isin(stats_filtered["sample"]),]
        logging.info("Successfully filtered samples from tables.\n")
    else:
        logging.info("No filtering requested, copying tables...")
        stats_filtered = stats
        metadata_filtered = metadata
    logging.info("Saving filtered tables...")
    stats_filtered.to_csv(snakemake.output.stats, index=False, header=True, sep = "\t")
    metadata_filtered.to_csv(snakemake.output.metadata, index=False)
    logging.info("Done!")
except Exception as e:
    log_file.write(f"Error: {e}\n")
    raise e