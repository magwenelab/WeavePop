import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

input_tsv=snakemake.input.tsv
output_bed=snakemake.output.bed

print("Reading GFF table...")
df = pd.read_csv(input_tsv, sep='\t', header=[0], dtype=str)

print("Converting to BED format...")
df = df[df['primary_tag'] == "gene"]
df = df[['accession', 'start', 'end', 'feature_id']]

print("Saving BED...")
df.to_csv(output_bed, sep='\t', index=False, header=False)
