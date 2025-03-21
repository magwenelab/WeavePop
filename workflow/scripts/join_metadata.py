import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

dataset_names = snakemake.params[0]
input=snakemake.input
output=snakemake.output

print("Getting dataset names...")
print("Reading metadata files and adding dataset names...")
dataframes = []
for file_path, string in zip(input, dataset_names):
    df = pd.read_csv(file_path, sep=",", header=0)
    df["dataset"] = string
    dataframes.append(df)
metadata_df= pd.concat(dataframes, ignore_index=True)
metadata_df = metadata_df.drop_duplicates()
metadata_df.to_csv(str(output), sep=",", index=False)
print("Done!")
