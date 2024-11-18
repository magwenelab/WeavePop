
import pandas as pd
from pathlib import Path
import logging

log_file=snakemake.log[0]
input=snakemake.input[0]
output=snakemake.output[0]

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading metadata file...")
    metadata = pd.read_csv(input, sep=",", header=0)
    lineages = list(set(metadata["lineage"]))
    if len(lineages) == 0:
        message = (
            "No lineages found. Exiting. "
        ) 
        raise ValueError(message)
    else:
        for lineage in lineages:
            path = Path(output, f"{lineage}.lineage")
            logging.info(f"Creating file: {path}..")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e

