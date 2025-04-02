import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

from pathlib import Path
import pandas as pd

input_cds = snakemake.input.cds
input_prots = snakemake.input.prots
output = snakemake.output.sequences

print("Create dictionary with lineage and paths for cds files...")
paths_cds = {}
for path in input_cds:
    lineage = Path(Path(path).stem).stem
    paths_cds[lineage] = path
    print("Available lineages and paths:")
    print(paths_cds)
    
print("Create dictionary with lineage and paths for protein files...")
paths_prots = {}
for path in input_prots:
    lineage = Path(Path(path).stem).stem
    paths_prots[lineage] = path
    print("Available lineages and paths:")
    print(paths_prots)
    
print("Cds files to concatenate:")
print(paths_cds)
print("Protein files to concatenate:")
print(paths_prots)


print("Reading sequences...")
cds = pd.concat([pd.read_csv(f, sep="\t") for f in paths_cds.values()])
prots = pd.concat([pd.read_csv(f, sep="\t") for f in paths_prots.values()])
print("Joining sequences...")
sequences = pd.concat([cds, prots])
print("Writing sequences...")
sequences.to_csv(output, sep="\t", index=False)
print("Done!")