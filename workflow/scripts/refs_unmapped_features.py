import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import os
import pandas as pd

genes_path=snakemake.input.genes

lineages_path=str(snakemake.params.refdir)
lineages = snakemake.params.lins

output_path=snakemake.output.tsv

print("Reading files...")
genes = pd.read_csv(genes_path, header=0, sep="\t")
print("Selecting columns of main reference...")
columns_to_select = ["seq_id", "start", "end", "primary_tag", "ID", "description" , "Name"]
existing_columns = [col for col in columns_to_select if col in genes.columns]
genes = genes[existing_columns]
print("Renaming columns...")
rename_dict = {
    "seq_id": "accession",
    "ID": "gene_id",
    "Name": "gene_name"
}
genes.rename(columns=rename_dict, inplace=True)

print("Reading list of unmapped features of each lineage and joining to annotation table...")
for lineage in lineages:
    unmapped_features_path = os.path.join(lineages_path, lineage, "unmapped_features.txt")
    unmapped_features = pd.read_csv(unmapped_features_path, sep="\t", header =None, names=["gene_id"])
    unmapped_features[lineage] = "unmapped"
    genes = genes.merge(unmapped_features, on="gene_id", how="left")

print("Filtering unmapped features...")
unmapped = genes.dropna(subset=lineages, how='all')

print("Saving table...")
unmapped.to_csv(output_path, sep="\t", index=False)
print("Done!")
