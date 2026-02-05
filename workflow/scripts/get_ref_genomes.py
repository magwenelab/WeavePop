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
ref_genomes = list(set(metadata["ref_genome"]))
if len(ref_genomes) == 0:
    message = (
        "No ref_genomes found. Exiting. "
    ) 
    raise ValueError(message)
else:
    for ref_genome in ref_genomes:
        path = Path(output, f"{ref_genome}.ref_genome")
        print(f"Creating file: {path}..")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    print("Done!")

