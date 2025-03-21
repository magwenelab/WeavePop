import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd
from pathlib import Path

input=snakemake.input[0]
output=snakemake.output[0]


print("Reading metadata file...")
metadata = pd.read_csv(input, sep=",", header=0)
lineages = list(set(metadata["lineage"]))
if len(lineages) == 0:
    message = (
        "No lineages found. Exiting. "
    ) 
    raise ValueError(message)
else:
    for lineage in lineages:
        path = Path(output, f"{lineage}.lineage")
        print(f"Creating file: {path}..")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    print("Done!")

