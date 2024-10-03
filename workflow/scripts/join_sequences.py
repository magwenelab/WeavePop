import pandas as pd
import logging

log_file = snakemake.log[0]
input_cds = snakemake.input.cds
input_prots = snakemake.input.prots
output = snakemake.output.sequences
# Configure logging
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading sequences...")
    cds = pd.concat([pd.read_csv(f, sep="\t") for f in input_cds])
    prots = pd.concat([pd.read_csv(f, sep="\t") for f in input_prots])
    logging.info("Joining sequences...")
    sequences = pd.concat([cds, prots])
    logging.info("Writing sequences...")
    sequences.to_csv(output, sep="\t", index=False)
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e