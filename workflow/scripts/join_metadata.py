import pandas as pd
import logging

log_file = snakemake.log[0]
dataset_names = snakemake.params[0]
input=snakemake.input
output=snakemake.output

logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')

try:
    logging.info("Getting dataset names...")
    logging.info("Reading metadata files and adding dataset names ...")
    dataframes = []
    for file_path, string in zip(input, dataset_names):
        df = pd.read_csv(file_path, sep=",", header=0)
        df["dataset"] = string
        dataframes.append(df)
    metadata_df= pd.concat(dataframes, ignore_index=True)
    metadata_df = metadata_df.drop_duplicates()
    metadata_df.to_csv(str(output), sep=",", index=False)
    logging.info("Done!")
except Exception as e:
    logging.error(f"Error: {e}")
    raise e
