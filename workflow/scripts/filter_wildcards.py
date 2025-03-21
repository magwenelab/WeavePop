import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd
from pathlib import Path

input=snakemake.input[0]
output_samples=snakemake.output[0]
output_lineages=snakemake.output[1]

print("Reading file...")
metadata = pd.read_csv(input, header=0)
sample_names = list(metadata["sample"])
lineage_names = set(metadata["lineage"])
if len(sample_names) == 0:
    message = (
        "The quality filter removed all samples! "
        "There is nothing to analyze. Exiting. "
        "See 04.Intermediate_files/02.Dataset/depth_quality/unfiltered_mapping_stats.tsv "
        "to check the quality warning of each sample."
    ) 
    raise ValueError(message)
else:
    print(f"Number of samples: {len(sample_names)}")
    for sample_name in sample_names:
        path = Path(output_samples,f"{sample_name}.txt")
        print(f"Creating file: {path}...")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    print(f"Number of lineages: {len(lineage_names)}")
    for lineage in lineage_names:
        path = Path(output_lineages,f"{lineage}.txt")
        print(f"Creating file: {path}...")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
print("Done!")


