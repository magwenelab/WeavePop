import duckdb
import pandas as pd
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from enum import Enum


def fasta_to_df(fasta, sample, seq_type):
    df_all = pd.DataFrame()
    for record in SeqIO.parse(fasta, "fasta"):
        transcript_id = record.id
        seq_description = record.description
        locus = seq_description.split(" ")[1]
        locus = locus.replace("gene", "locus")
        accession = seq_description.split(" ")[2]
        accession = accession.replace("seq_id", "accession")
        new_description = sample + "|" + transcript_id + " " + locus + " " + accession
        seq = str(record.seq)
        df = pd.DataFrame({"sample":sample, "transcript_id":transcript_id, "seq":seq, "seq_type":seq_type,  "seq_description":new_description}, index=[0])
        df_all = pd.concat([df_all, df])
    return df_all

def sample_type_not_in_db(con, sample, seq_type):
    result = con.execute(f"SELECT * FROM sequences WHERE sample = '{sample}' AND seq_type = '{seq_type}'").fetch_df()
    if result.shape[0] > 0:
        return False
    else:
        return True

def populate_fasta_db(db, fasta, sample, seq_type):
    print("Creating dataframe")
    seq_type = seq_type.upper()
    df = fasta_to_df(fasta, sample, seq_type)   
    con = duckdb.connect(db)
    result = con.execute("SHOW TABLES").fetch_df()
    if "sequences" not in result["name"].to_list():
        print("Creating table and populating it with dataframe")
        con.register("df_sequences", df)
        con.execute("CREATE TABLE sequences AS SELECT * FROM df_sequences")
        con.close()
    elif sample_type_not_in_db(con, sample, seq_type):
        print("Populating existing table with dataframe")
        con.register("df_sequences", df)
        con.execute("INSERT INTO sequences SELECT * FROM df_sequences")
        con.close()
    else:
        print(f"Sample {sample} with sequence type {seq_type} already in database")
        con.close()
    print("Done")

@click.command()
@click.option("--db", "-d", required=True, help="Path to the database file")
@click.option("--fasta", "-f", required=True, help="Path to the fasta file")
@click.option("--sample", "-s", required=True, help="Sample name")
@click.option("--seq_type", "-t", required=True, help="Sequence type")
def main(db, fasta, sample, seq_type):
    populate_fasta_db(db, fasta, sample, seq_type)
    
if __name__ == "__main__":
    main()        