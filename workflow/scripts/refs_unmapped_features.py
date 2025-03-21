import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd
import os

metadata_path=snakemake.input[0]
lineages_path=snakemake.params[0]
genes_path=snakemake.input[1]
output_path=snakemake.output[0]

print("Reading files...")
metadata = pd.read_csv(metadata_path, header=0)
genes = pd.read_csv(genes_path, header=0, sep="\t")
print("Getting lineages...")
lineages = metadata["lineage"].unique()
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
