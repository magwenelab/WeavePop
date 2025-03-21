import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

input=snakemake.input
output=snakemake.output[0]

print("Reading and concatenating files...")
file = pd.concat([pd.read_csv(f, sep="\t") for f in input])
file.to_csv(output, sep="\t", index=False)
print("Done!")

