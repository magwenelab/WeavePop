import pandas as pd
import duckdb
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

cwd = os.getcwd()
def df_to_seqrecord(df):
    records = []
    for index, row in df.iterrows():
        seq = Seq(row['seq'])
        record = SeqRecord(seq, id=f"{row['sample']}|{row['transcript_id']}", description=row['seq_description'])
        records.append(record)
    return records

def get_sequences_of_gene(gene_name,db='database/desjardins.db', seqtype='DNA'):
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    query = f"""
        SELECT *
        FROM sequences
        WHERE transcript_id IN (
            SELECT DISTINCT feature_id
            FROM gff
            WHERE parent = '{gene_name}'
        )
        AND seq_type = '{seqtype}'"""
    result = con.execute(query).fetchdf()
    con.close()
    if seqtype == 'PROTEIN':
        ext = 'faa'
    elif seqtype == 'DNA':
        ext = 'fna'
    file_name = f'{gene_name}.{ext}'
    print(f'Saving result to fasta file {file_name}')
    SeqIO.write(df_to_seqrecord(result), file_name, 'fasta')
    
def get_samples_with_variants_in_gene(gene_name,db='database/desjardins.db', impact=None):
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    query = f"""
        SELECT DISTINCT sample
        FROM presence 
        WHERE var_id IN (
            SELECT var_id 
            FROM effects 
            WHERE locus = '{gene_name}'"""
    if impact:
        query += f" AND impact = '{impact}'"
    query += ")"
    result = con.execute(query).fetchdf()
    con.close()
    return result

def get_sequences_of_gene_with_variants(gene_name,db='database/desjardins.db',seqtype='DNA', impact=None):
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    samples = get_samples_with_variants_in_gene(gene_name,db, impact)
    samples_tuple = tuple(samples['sample'])
    
    if not samples_tuple:
        print("There are no samples that meet your specifications")
        return
    
    query = f"""
        SELECT *
        FROM sequences
        WHERE sample IN {samples_tuple}
            AND seq_type = '{seqtype}'
            AND transcript_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE parent = '{gene_name}')"""
                
    result = con.execute(query).fetchdf()
    con.close()
    if seqtype == 'PROTEIN':
        ext = 'faa'
    elif seqtype == 'DNA':
        ext = 'fna'
        
    if impact:
        file_name = f'{gene_name}_{impact}.{ext}'
    else:
        file_name = f'{gene_name}.{ext}'
    print(f'Saving result to fasta file {file_name}')
    SeqIO.write(df_to_seqrecord(result), file_name, 'fasta')

def get_effects_of_variants_in_gene(gene_name,db='database/desjardins.db',impact=None):
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    query = f"""
        SELECT *
        FROM effects 
        WHERE locus = '{gene_name}'"""
    if impact:
        query += f" AND impact = '{impact}'"
    result = con.execute(query).fetchdf()
    con.close()
    return result

def get_variants_with_effect_in_gene(gene_name,db='database/desjardins.db',impact=None):
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    vars = get_effects_of_variants_in_gene(gene_name,db,impact)
    vars_tuple = tuple(vars['var_id'])
    query = f"""
        SELECT *
        FROM presence 
        WHERE var_id IN {vars_tuple}""" 
    result = con.execute(query).fetchdf()
    con.close()
    return result

def count_variants_in_gene(gene_name,db='database/desjardins.db',impact=None):
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    vars = get_effects_of_variants_in_gene(gene_name,db,impact)
    vars_tuple = tuple(vars['var_id'])
    query = f"""
        SELECT sample, lineage, COUNT(var_id) num_vars
        FROM presence 
        WHERE var_id IN {vars_tuple}
        GROUP BY sample, lineage""" 
    result = con.execute(query).fetchdf()
    con.close()
    return result

def get_lineage_of_samples_with_variants_in_gene(gene_name,db='database/desjardins.db',impact=None):
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    vars = get_effects_of_variants_in_gene(gene_name,db,impact)
    lins_tuple = tuple(vars['lineage'])
    query = f"""
        SELECT sample, lineage
        FROM samples 
        WHERE lineage IN {lins_tuple}""" 
    result = con.execute(query).fetchdf()
    con.close()
    return result

def get_summary_variants_in_gene(gene_name, db='database/desjardins.db', impact=None):
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    counts = count_variants_in_gene(gene_name,db,impact)
    vars = get_effects_of_variants_in_gene(gene_name,db,impact)
    samples_lins = get_lineage_of_samples_with_variants_in_gene(gene_name,db,impact)
    con.close()
    summary = counts.groupby('lineage').agg({'num_vars': 'mean'}).reset_index()
    summary.rename(columns={'num_vars': 'mean_vars_per_sample'}, inplace=True)
    summary['num_vars_in_lineage'] = vars.groupby('lineage')['var_id'].nunique().values
    summary['num_samples_with_vars'] = counts.groupby('lineage').size().values
    summary['total_samples_in_lineage'] = samples_lins.groupby('lineage').size().values
    return summary
