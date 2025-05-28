
import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

input_stats=snakemake.input[0]
input_chromosomes=snakemake.input[1]
filter=snakemake.params.filter
output_metadata=snakemake.output.metadata
output_chromosomes=snakemake.output.chromosomes
metadata = snakemake.params.metadata


print("Reading tables...")
stats = pd.read_csv(input_stats, sep="\t", header=0)
chromosomes = pd.read_csv(input_chromosomes, sep= ",", header=0)
if filter:
    print("Filtering samples...")
    stats_filtered = stats[stats["quality_warning"].isna()]
    metadata_filtered = metadata.loc[metadata["sample"].isin(stats_filtered["sample"]),]
    chromosomes_filtered = chromosomes.loc[chromosomes["lineage"].isin(metadata_filtered["lineage"]),]
    print("Successfully filtered samples from tables.\n")
else:
    print("No filtering requested, copying tables...")
    metadata_filtered = metadata
    chromosomes_filtered = chromosomes

print("Cleaning column names...")
metadata_filtered.columns = metadata_filtered.columns.str.lower()
metadata_filtered.columns = metadata_filtered.columns.str.replace(' ', '_')
print("Adding dataset column...")
metadata_filtered.loc[:, 'dataset'] = "X"
    
print("Saving filtered tables...")
metadata_filtered.to_csv(output_metadata, index=False)
chromosomes_filtered.to_csv(output_chromosomes, index=False)
print("Done!")