import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

input_cds = snakemake.input.cds
input_prots = snakemake.input.prots
output = snakemake.output.sequences

print("Reading sequences...")
cds = pd.concat([pd.read_csv(f, sep="\t") for f in input_cds])
prots = pd.concat([pd.read_csv(f, sep="\t") for f in input_prots])
print("Joining sequences...")
sequences = pd.concat([cds, prots])
print("Writing sequences...")
sequences.to_csv(output, sep="\t", index=False)
print("Done!")