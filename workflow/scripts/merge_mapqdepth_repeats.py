import pandas as pd
import logging

log_file=snakemake.log[0]
input_depth=snakemake.input.depth
input_repeats=snakemake.input.repeats
output=snakemake.output

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading tables...")
    depth = pd.read_csv(input_depth, sep="\t", header=0)
    repeats = pd.read_csv(input_repeats, sep="\t", header=0)
    logging.info("Merging tables...")
    merged = pd.merge(depth, repeats, on="ID", how="left")
    logging.info("Saving merged table...")
    merged.to_csv(output[0], index=False, header=True, sep = "\t")
    logging.info("Done!")
except Exception as e:
    log_file.write(f"Error: {e}\n")
    raise e