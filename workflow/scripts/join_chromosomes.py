import pandas as pd
import logging

log_file = snakemake.log[0]

# Configure logging
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading and concatenting chromosomes...")
    chromosomes_df = pd.concat([pd.read_csv(f, sep=",", header=0) for f in snakemake.input])
    chromosomes_df = chromosomes_df.drop_duplicates()
    logging.info("Saving chromosomes file...")
    chromosomes_df.to_csv(str(snakemake.output), sep=",", index=False)
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e