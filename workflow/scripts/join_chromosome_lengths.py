import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

chrom_names=snakemake.input.chrom_names
chrom_lengths=snakemake.input.chrom_lengths 
output=snakemake.output

print("Reading and concatenting chromosome lengths...")
chrom_lengths_df = pd.concat([pd.read_csv(f, sep="\t", names = ["lineage", "accession", "length"]) for f in chrom_lengths])
print("Reading chromosome names...")
chrom_names_df = pd.read_csv(chrom_names, sep=",", header=0)
print("Merging chromosome names and lengths...")
chromosomes_df =  chrom_names_df.merge(chrom_lengths_df, how="left", left_on=["accession", "lineage"], right_on=["accession", "lineage"])
print("Saving chromosomes file...")
chromosomes_df.to_csv(str(output), sep="\t", index=False)
print("Done!")