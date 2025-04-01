import pandas as pd
import click
from Bio import SeqIO

@click.command()
@click.option("--fasta", "-f", required=True, type=click.Path(), help="Path to the fasta file")
@click.option("--sample", "-s", required =False, type=str, help="Sample name")
@click.option("--lineage", "-l", required =False, type=str, help="Lineage name")
@click.option("--seq_type", "-t", required=True, type=str, help="Sequence type")
@click.option("--output", "-o", required=True, type=click.Path(), help="Path to the output file")
def fasta_to_csv(fasta, sample,lineage, seq_type, output):
    df_all = pd.DataFrame()
    for record in SeqIO.parse(fasta, "fasta"):
        print("Extracting information for " + record.id + "...")
        transcript_id = record.id
        seq_description = record.description
        gene_id = seq_description.split(" ")[1]
        gene_id = gene_id.replace("gene", "gene_id")
        accession = seq_description.split(" ")[2]
        accession = accession.replace("seq_id", "accession")
        seq = str(record.seq)
        print("Joining information into a dataframe...")
        if sample is not None:
            new_description = sample + "|" + transcript_id + " " + gene_id + " " + accession
            df = pd.DataFrame({"sample":sample, "transcript_id":transcript_id, "seq":seq, "seq_type":seq_type,  "seq_description":new_description}, index=[0])
        elif lineage is not None:
            new_description = lineage + "|" + transcript_id + " " + gene_id + " " + accession
            df = pd.DataFrame({"lineage":lineage, "transcript_id":transcript_id, "seq":seq, "seq_type":seq_type,  "seq_description":new_description}, index=[0])
        if df.isnull().values.any() or (df == '').any().any():
            raise ValueError("Empty or NaN value found in data for sample " + sample + " and transcript " + transcript_id + ". Check fasta file.")
        print("Adding record to the dataframe with all records...")
        df_all = pd.concat([df_all, df])
    print("Saving dataframe...")
    df_all.to_csv(output, index=False, header=True)

if __name__ == "__main__":
    fasta_to_csv()    