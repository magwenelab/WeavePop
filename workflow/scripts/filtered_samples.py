
import pandas as pd
from pathlib import Path
import logging

log_file=snakemake.log[0]
input=snakemake.input[0]
output=snakemake.output[0]

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading file...")
    metadata = pd.read_csv(input, header=0)
    sample_names = list(metadata["sample"])
    for sample_name in sample_names:
        path = Path(output,f"{sample_name}.txt")
        logging.info(f"Creating file: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e

