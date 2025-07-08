import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

chrom_names=snakemake.input[0]
chrom_lengths=snakemake.input[1]
output=snakemake.output

print("Reading and concatenting chromosome lengths...")
chrom_lengths_df = pd.read_csv(chrom_lengths, sep="\t", header=None, names = ["accession", "length", "lineage"])
print("Reading chromosome names...")
chrom_names_df = pd.read_csv(chrom_names, sep=",", header=0)
print("Merging chromosome names and lengths...")
chromosomes_df =  chrom_lengths_df.merge(chrom_names_df, how="left", on=["accession", "lineage"])
chromosomes_df = chromosomes_df[["accession", "chromosome", "length", "lineage"]]
print("Saving chromosomes file...")
chromosomes_df.to_csv(str(output), sep=",", index=False)
print("Done!")