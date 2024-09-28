import pandas as pd
import logging

log_file = snakemake.log[0]

# Configure logging
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading and concatenating files...")
    file = pd.concat([pd.read_csv(f, sep="\t") for f in snakemake.input])
    file.to_csv(snakemake.output[0], sep="\t", index=False)
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e

