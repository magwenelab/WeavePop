
import pandas as pd
from pathlib import Path
import logging

log_file = snakemake.log[0]

# Configure logging
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading file...")
    metadata = pd.read_csv(snakemake.input[0], header=0)
    lineages = list(metadata["lineage"])
    for lineage in lineages:
        path = Path(snakemake.output[0],f"{lineage}.txt")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e

