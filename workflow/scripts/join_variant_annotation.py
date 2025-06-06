import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

from pathlib import Path
import pandas as pd

input_effects=snakemake.input.effects
input_variants=snakemake.input.variants
input_lofs=snakemake.input.lofs
input_nmds=snakemake.input.nmds
input_presence=snakemake.input.presence
output_effects=snakemake.output.effects
output_variants=snakemake.output.variants
output_lofs=snakemake.output.lofs
output_nmds=snakemake.output.nmds
output_presence=snakemake.output.presence

print("Reading and concatenating variant annotation files...")
effects = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in input_effects])
variants = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in input_variants])
lofs = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in input_lofs])
nmds = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in input_nmds])
presence = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in input_presence])
print("Saving variant annotation files...")
effects.to_csv(output_effects, sep = "\t", index=False)
variants.to_csv(output_variants, sep="\t", index=False)
lofs.to_csv(output_lofs, sep="\t", index=False)
nmds.to_csv(output_nmds, sep="\t", index=False)
presence.to_csv(output_presence, sep="\t", index=False)
print("Done!")