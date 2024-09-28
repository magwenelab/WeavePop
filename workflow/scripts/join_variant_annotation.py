from pathlib import Path
import pandas as pd
import logging

log_file = snakemake.log[0]

# Configure logging
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading and concatenating variant annotation files...")
    effects = pd.concat([pd.read_csv(f, sep="\t") for f in snakemake.input.effects])
    variants = pd.concat([pd.read_csv(f, sep="\t") for f in snakemake.input.variants])
    lofs = pd.concat([pd.read_csv(f, sep="\t") for f in snakemake.input.lofs])
    nmds = pd.concat([pd.read_csv(f, sep="\t") for f in snakemake.input.nmds])
    presence = pd.concat([pd.read_csv(f, sep="\t") for f in snakemake.input.presence])
    logging.info("Saving variant annotation files...")
    effects.to_csv(snakemake.output.effects, sep = "\t", index=False)
    variants.to_csv(snakemake.output.variants, sep="\t", index=False)
    lofs.to_csv(snakemake.output.lofs, sep="\t", index=False)
    nmds.to_csv(snakemake.output.nmds, sep="\t", index=False)
    presence.to_csv(snakemake.output.presence, sep="\t", index=False)
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e