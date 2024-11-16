import pandas as pd
import duckdb
import click
import os
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

cwd = os.getcwd()

def sequences(db, dataset=None, seq_type='DNA', sample=None, strain=None, lineage=None, gene_id=None, gene_name=None):
    if gene_name and gene_id:
        raise ValueError("Only one of Gene names or Gene IDs should be provided.")
    else:
        pass
    if sample and strain or sample and lineage or strain and lineage:
        raise ValueError("Only one of Sample IDs, Strains or Lineage should be provided.")
    else:
        pass
    
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    query_seq_type = f"""
        SELECT DISTINCT seq_type
        FROM sequences
        """
    seq_type_df = con.execute(query_seq_type).fetchdf()
    seq_type_list = tuple(seq_type_df['seq_type'].tolist())
    if seq_type not in seq_type_list:
        print(f"seq_type must be one of {seq_type_list}. Exiting.")
        sys.exit()
    
    query_dataset = f"""
        SELECT DISTINCT dataset
        FROM samples
        """
    dataset_df = con.execute(query_dataset).fetchdf()
    dataset_list = tuple(dataset_df['dataset'].tolist())
    
    if dataset is None:
        dataset = dataset_list
    else:
        dataset_input = tuple(dataset.split(','))
        if all(d in dataset_list for d in dataset_input):
            dataset = dataset_input
        else:
            if any(d not in dataset_list for d in dataset_input):
                dataset_absent = tuple(d for d in dataset_input if d not in dataset_list)
                print(f"{dataset_absent} not found in the database. Excluding.")
                if len(dataset_absent) == len(dataset_input):
                    print("No valid dataset found. Exiting.")
                    sys.exit()
                else:
                    dataset = tuple(d for d in dataset_input if d in dataset_list)
        

    if gene_name is not None:
        query_gene_name = f"""
            SELECT DISTINCT gene_name
            FROM gff
            """
        gene_name_df = con.execute(query_gene_name).fetchdf()
        gene_name_list = tuple(gene_name_df['gene_name'].tolist())
        gene_name_input = tuple(gene_name.split(','))
        if all(d in gene_name_list for d in gene_name_input):
            gene_name = gene_name_input
        else:
            if any(d not in gene_name_list for d in gene_name_input):
                gene_name_absent = tuple(d for d in gene_name_input if d not in gene_name_list)
                print(f"{gene_name_absent} not found in the database. Excluding.")
                if len(gene_name_absent) == len(gene_name_input):
                    print("No valid gene name found. Exiting.")
                    sys.exit()
                else:
                    gene_name = tuple(d for d in gene_name_input if d in gene_name_list)
        query_gene_id = f"""
            SELECT *
            FROM gff
            WHERE gene_name IN {gene_name}
            """
        gene_id_df = con.execute(query_gene_id).fetchdf()
        gene_id = tuple(gene_id_df['gene_id'].unique().tolist())
    elif gene_id is not None:
        query_gene_id = f"""
            SELECT DISTINCT gene_id
            FROM gff
            """
        gene_id_df = con.execute(query_gene_id).fetchdf()
        gene_id_list = tuple(gene_id_df['gene_id'].tolist())     
        gene_id_input = tuple(gene_id.split(','))
        if all(d in gene_id_list for d in gene_id_input):
            gene_id = gene_id_input
        else:
            if any(d not in gene_id_list for d in gene_id_input):
                gene_id_absent = tuple(d for d in gene_id_input if d not in gene_id_list)
                print(f"{gene_id_absent} not found in the database. Excluding.")
                if len(gene_id_absent) == len(gene_id_input):
                    print("No valid gene ID found. Exiting.")
                    sys.exit()
                else:
                    gene_id = tuple(d for d in gene_id_input if d in gene_id_list)
        
    if strain is not None:
        query_strain = f"""
            SELECT DISTINCT strain
            FROM samples
            """
        strain_df = con.execute(query_strain).fetchdf()
        strain_list = tuple(strain_df['strain'].tolist())
        strain_input = tuple(strain.split(','))
        if all(d in strain_list for d in strain_input):
            strain = strain_input
        else:
            if any(d not in strain_list for d in strain_input):
                strain_absent = tuple(d for d in strain_input if d not in strain_list)
                print(f"{strain_absent} not found in the database. Excluding.")
                if len(strain_absent) == len(strain_input):
                    print("No valid strain found. Exiting.")
                    sys.exit()
                else:
                    strain = tuple(d for d in strain_input if d in strain_list)
        query_sample = f"""
            SELECT sample
            FROM samples
            WHERE strain IN {strain}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'].tolist())
    elif lineage is not None:
        query_lineage = f"""
            SELECT DISTINCT lineage
            FROM samples
            """
        lineage_df = con.execute(query_lineage).fetchdf()
        lineage_list = tuple(lineage_df['lineage'].tolist())
        lineage_input = tuple(lineage.split(','))
        if all(d in lineage_list for d in lineage_input):
            lineage = lineage_input
        else:
            if any(d not in lineage_list for d in lineage_input):
                lineage_absent = tuple(d for d in lineage_input if d not in lineage_list)
                print(f"{lineage_absent} not found in the database. Excluding.")
                if len(lineage_absent) == len(lineage_input):
                    print("No valid lineage found. Exiting.")
                    sys.exit()
                else:
                    lineage = tuple(d for d in lineage_input if d in lineage_list)
        query_sample = f"""
            SELECT sample
            FROM samples
            WHERE lineage IN {lineage}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'].tolist())
    elif sample is not None:
        sample_input = tuple(sample.split(','))
        query_sample = f"""
            SELECT DISTINCT sample
            FROM samples
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample_list = tuple(sample_df['sample'].tolist())
        if all(d in sample_list for d in sample_input):
            sample = sample_input
        else:
            if any(d not in sample_list for d in sample_input):
                sample_absent = tuple(d for d in sample_input if d not in sample_list)
                print(f"{sample_absent} not found in the database. Excluding.")
                if len(sample_absent) == len(sample_input):
                    print("No valid sample found. Exiting.")
                    sys.exit()
                else:
                    sample = tuple(d for d in sample_input if d in sample_list)
   
    query = f"""
            SELECT samples.dataset, samples.strain, samples.lineage, 
                sequences.sample, sequences.transcript_id, sequences.seq, 
                chromosome_names.chromosome, chromosome_names.accession,
                gff.gene_name, gff.gene_id
            FROM sequences
            JOIN samples ON sequences.sample = samples.sample
            JOIN gff ON sequences.transcript_id = gff.feature_id AND samples.lineage = gff.lineage
            JOIN chromosome_names ON gff.accession = chromosome_names.accession
            WHERE samples.dataset IN {dataset}"""
    if gene_id and not sample:
        query += f"""
            AND transcript_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )
            AND seq_type = '{seq_type}'"""
    elif gene_id and sample:
        query += f"""
            AND transcript_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )
            AND sequences.sample IN {sample}
            AND seq_type = '{seq_type}'"""
    elif sample:
        query += f"""
            AND sequences.sample IN {sample}
            AND seq_type = '{seq_type}'"""
    result = con.execute(query).fetchdf()
    con.close()
    return result

def df_to_seqrecord(df):
    records = []
    for index, row in df.iterrows():
        seq = Seq(row['seq'])
        record = SeqRecord(seq, id=f"{row['strain']}|{row['transcript_id']}", description=f"sample={row['sample']} gene_id={row['gene_id']} gene_name={row['gene_name']} chromosome={row['chromosome']} accession={row['accession']}")
        records.append(record)
    return records

def seqrecord_to_text(records):
    output = StringIO()
    SeqIO.write(records, output, "fasta")
    fasta_text = output.getvalue()
    output.close()
    return fasta_text

@click.command()
@click.option('--db', help='Path to the database file')
@click.option('--gene_id', default=None, help='Comma separated list of gene IDs')
@click.option('--gene_name', default=None, help='Comma separated list of gene names')
@click.option('--dataset', default=None, help='Comma separated list of dataset names')
@click.option('--sample', default=None, help='Comma separated list of sample IDs')
@click.option('--strain', default=None, help='Comma separated list of strain names')
@click.option('--lineage', default=None, help='Comma separated list of lineage names')
@click.option('--seq_type', default='DNA', help='Sequence type. Options are DNA and PROTEIN', show_default=True)
@click.option('--filename', default=None, help='Output filename. Printed to standard output if not provided.')
def main(db, gene_id, gene_name, dataset, sample, strain, lineage, seq_type, filename):
    result = sequences(db, gene_id=gene_id, gene_name=gene_name, dataset=dataset, sample=sample, strain=strain, lineage=lineage, seq_type=seq_type)
    records = df_to_seqrecord(result)
    fasta_text = seqrecord_to_text(records)
    if filename is None:
        print(fasta_text)
    else:
        with open(filename, 'w') as f:
            f.write(fasta_text)

if __name__ == '__main__':
    main()