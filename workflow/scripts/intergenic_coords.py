import pandas as pd
import logging

log_file=snakemake.log[0]
input_seqs=snakemake.input[0]
input_gff=snakemake.input[1]
output=snakemake.output[0]

# Configure logging
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Reading gff tsv...")
    gff_df = pd.read_csv(input_gff, sep="\t", header=0)
    logging.info("Selecting and renaming columns...")
    gff_df = gff_df[['ID', 'seq_id', 'start', 'end']]
    gff_df = gff_df.rename(columns={"ID": "transcript_id"})
    logging.info("Reading intergenic csv...")
    intergenic_df = pd.read_csv(input_seqs, sep=",", header=0)
    logging.info("Merging dataframes...")
    merged = pd.merge(intergenic_df, gff_df, on="transcript_id", how="left")
    logging.info("Saving merged file...")
    merged.to_csv(output, sep=",", index=False)
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e