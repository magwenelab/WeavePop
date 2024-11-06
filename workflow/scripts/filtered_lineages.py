
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
    lineages = list(metadata["lineage"])
    if len(lineages) == 0:
        message = (
            "The quality filter removed all samples! "
            "There is nothing to analyze. Exiting. "
            "See 04.Intermediate_files/02.Dataset/depth_quality/unfiltered_mapping_stats.tsv "
            "to check the quality warning of each sample."
        ) 
        raise ValueError(message)
    else:
        for lineage in lineages:
            path = Path(output,f"{lineage}.txt")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e

