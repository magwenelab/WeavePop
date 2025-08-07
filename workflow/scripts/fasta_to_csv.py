import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import pandas as pd
from Bio import SeqIO
    
fasta=snakemake.input.fa
seq_type=snakemake.params.seq_type

if snakemake.wildcards.get('sample') is not None:
    sample = snakemake.wildcards.sample
else:
    sample = None
if snakemake.wildcards.get('lineage') is not None:
    lineage = snakemake.wildcards.lineage
else:
    lineage = None

output=snakemake.output.csv

df_all = pd.DataFrame()
for record in SeqIO.parse(fasta, "fasta"):
    print("Extracting information for " + record.id + "...")
    feature_id = record.id
    seq_description = record.description
    gene_id = seq_description.split(" ")[1]
    gene_id = gene_id.replace("gene", "gene_id")
    accession = seq_description.split(" ")[2]
    accession = accession.replace("seq_id", "accession")
    seq = str(record.seq)
    print("Joining information into a dataframe...")
    if sample is not None:
        new_description = sample + "|" + feature_id + " " + gene_id + " " + accession
        df = pd.DataFrame({"sample":sample, "feature_id":feature_id, "seq_type":seq_type, "seq":seq, "seq_description":new_description}, index=[0])
    elif lineage is not None:
        new_description = lineage + "|" + feature_id + " " + gene_id + " " + accession
        df = pd.DataFrame({"lineage":lineage, "feature_id":feature_id, "seq_type":seq_type, "seq":seq, "seq_description":new_description}, index=[0])
    if df.isnull().values.any() or (df == '').any().any():
        raise ValueError("Empty or NaN value found in data for sample " + sample + " and transcript " + feature_id + ". Check fasta file.")
    print("Adding record to the dataframe with all records...")
    df_all = pd.concat([df_all, df])
print("Saving dataframe...")
df_all.to_csv(output, index=False, header=True)