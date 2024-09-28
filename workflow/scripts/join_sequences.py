import pandas as pd
import logging

log_file = snakemake.log[0]

# Configure logging
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading sequences...")
    cds = pd.concat([pd.read_csv(f, sep="\t") for f in snakemake.input.cds])
    prots = pd.concat([pd.read_csv(f, sep="\t") for f in snakemake.input.prots])
    logging.info("Joining sequences...")
    sequences = pd.concat([cds, prots])
    logging.info("Writing sequences...")
    sequences.to_csv(snakemake.output.sequences, sep="\t", index=False)
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e