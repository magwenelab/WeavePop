import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

input=snakemake.input
output=snakemake.output

print("Reading and concatenting chromosomes...")
chromosomes_df = pd.concat([pd.read_csv(f, sep=",", header=0) for f in input])
chromosomes_df = chromosomes_df.drop_duplicates()
print("Saving chromosomes file...")
chromosomes_df.to_csv(str(output), sep=",", index=False)
print("Done!")