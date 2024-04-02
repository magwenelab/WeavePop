#!/usr/bin/env python

from enum import Enum
import sqlite3
from pathlib import Path

import click

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


SCHEMA = """
CREATE TABLE IF NOT EXISTS samples (
    sample_id TEXT PRIMARY KEY,   -- systematic name (e.g. NCBI biosample identifier) 
    common_name TEXT,             -- common name (e.g. H99) 
    noop DEFAULT 0                -- convenient for upsert clauses where still want to return something
);
CREATE TABLE IF NOT EXISTS sequences (
    numid INTEGER PRIMARY KEY,  
    sample_id TEXT references samples,
    seq_id TEXT NOT NULL,
    seq TEXT NOT NULL,
    seq_type TEXT,
    seq_name TEXT,
    seq_description TEXT,
    seq_dbxrefs TEXT
);
"""


class SeqType(Enum):
    DNA = 1
    RNA = 2
    PROTEIN = 3


str2SeqType = {t.name:t for t in SeqType}



def create_tables(con, schema=SCHEMA):
    con.executescript(schema)
    con.commit()

def add_sample(con, sample_id, common_name=None):
    sid = con.execute("""insert into samples(sample_id, common_name) values(?, ?) 
                      on conflict(sample_id) do update set noop = 1 returning sample_id""", 
                      [sample_id, common_name]).fetchone()[0]    
    con.commit()
    return sid

def populate_fasta_db(outfile, infasta, sample_id, common_name=None, seq_type = None, timeout=60):
    con = sqlite3.connect(outfile,timeout=timeout)
    create_tables(con)
    sample_id = add_sample(con, sample_id, common_name)
    for record in SeqIO.parse(infasta, "fasta"):
        seq_id = record.id
        seq_name = record.name
        seq_description = record.description
        seq_dbxrefs = ";".join(record.dbxrefs)
        match seq_type:
            case None:
                stype = seq_type
            case _:
                stype = seq_type.name
        seq = str(record.seq)

        con.execute("""INSERT INTO 
                    sequences(numid, sample_id, seq_id, seq, seq_type, 
                              seq_name, seq_description, seq_dbxrefs)
                    VALUES(?, ?, ?, ?, ?, ?, ?, ?)""",
                    [None, sample_id, seq_id, seq, stype, 
                     seq_name, seq_description, seq_dbxrefs])
    con.commit()
    return con


def dict_factory(cursor, row):
    fields = [column[0] for column in cursor.description]
    return {key: value for key, value in zip(fields, row)}

def query_results(con, query):
    con.row_factory = dict_factory
    result = con.execute(query)
    while True:
        nextrow = result.fetchone()
        if not nextrow:
            break
        yield nextrow

def toSeqRecord(d, use_sample=True):
    if use_sample:
        id = f"{d['sample_id']}|{d['seq_id']}"
    else:
        id = f"{d['seq_id']}"
    trimmed_description = d["seq_description"].split(maxsplit=1)[-1]
    return SeqRecord(Seq(d["seq"]),
                     id = id,
                     name = d["seq_name"],
                     description = f"sample={d['sample_id']} {trimmed_description}",
                     dbxrefs = d["seq_dbxrefs"].split(';'))
    
    


@click.group()
def cli():
    """
    Tool/library for working with a Sqlite database with the following schema:

    \b
    samples (
        sample_id TEXT PRIMARY KEY,   -- systematic name (e.g. NCBI biosample identifier) 
        common_name TEXT,             -- common name (e.g. H99) 
        noop DEFAULT 0                -- convenient for upsert clauses where still want to return something
    );
    sequences (
        numid INTEGER PRIMARY KEY,  
        sample_id TEXT references samples,
        seq_id TEXT NOT NULL,
        seq TEXT NOT NULL,
        seq_type TEXT,
        seq_name TEXT,
        seq_description TEXT,
        seq_dbxrefs TEXT
    );
    """
    pass


@cli.command()
@click.option("--common_name", "-c",
              type=str,
              help="Common name associated with systematic sample name.")
@click.option("--seqtype", "-t",
              type=click.Choice(['DNA', 'RNA', 'PROTEIN', 'None'], case_sensitive=False),
              help="Indicates type of sequence information present in input FASTA file (DNA, RNA, PROTEIN, None (NULL))")
@click.option("--timeout",
              type=int,
              default=60,
              help="Timeout time (seconds) the Sqlite database should wait "
              " when a table is locked before raising an error.")
@click.argument("dbfile", type=click.Path(dir_okay=False,exists=False))
@click.argument("infasta", type=click.Path(dir_okay=False,exists=True))
@click.argument("sample_id", type=str)
def populate_db(dbfile, infasta, sample_id, common_name, seqtype=None, timeout=60):
    """
    Populate a database (DBFILE) with  records from the input FASTA file (INFASTA), 
    tagging them as associated with the given sample (SAMPLE_ID).
    """
    if seqtype:
        seqtype = str2SeqType[seqtype.upper()]
    con = populate_fasta_db(dbfile, infasta, sample_id, common_name, seqtype, timeout=timeout)
    con.close()


@cli.command()
@click.argument("db", type=click.Path(dir_okay=False, exists=True))
@click.argument("output", type=click.File('w'))
@click.option("--query", type=str)
@click.option("--query_file", type=click.Path(exists=True))
def run_query(db, output, query="", query_file=None):
    """ Run a given SQL query against the Sqlite database. 

    \b
    If --query is specified, the query is given on the command line.
    If --query_file is specified, the query is read from the file.
    If both --query and --query_file are specified, --query_file takes precedence.
    """
    if query_file:
        with open(query_file) as f:
            query = f.read()

    con = sqlite3.connect(db)
    for result in query_results(con, query):
        SeqIO.write(toSeqRecord(result), output, "fasta")
    con.close()



@cli.command()
@click.option("--seqid",
              multiple=True,
              help="Sequence ID to search on. If not specified will return all sequence IDs. Accepts glob patterns.")
@click.option("--sample",
              multiple=True,
              help="Sample ID to search on. Can be specified multiple times. "
                "If not specified returns matches for all samples.")
@click.option("--samplefile",
              type = click.Path(exists=True),
              help="Text file giving names of samples to include in search (one sample per line).")
@click.option("--seqtype", "-t",
              type=click.Choice(['DNA', 'RNA', 'PROTEIN',  'None'], case_sensitive=False),
              default="None",
              help="Indicates type of sequence information present in input FASTA file (DNA, RNA, PROTEIN, None (NULL))")
@click.option("--showquery", 
              is_flag=True,
              default=False,
              help="This option shows the query that would be constructed but doesn't run the query.")
@click.argument("db", type=click.Path(dir_okay=False, exists=True))
@click.argument("output", type=click.File('w'))
def lookup(db, output, sample, seqid, seqtype, samplefile=None, showquery=False):
    """Query the given sample-sequence Sqlite database returning sequences for a seq_id of interest.

    If --sample or --samplefile not specified, returns all samples with given seq_id of interest. 
    --sample can specified one or multiple times to manually specify one or a few samples. Alternately
    a file can be specified with --samplefile to give a longer list of samples to include in the search.
    """
    if samplefile:
        fromfile = []
        with open(samplefile) as f:
            for line in f:
                fromfile.append(line.strip())
        sample = list(sample.strip()) + fromfile
    if not sample:
        samplecond = "(SELECT sample_id FROM samples)"
    else:
        samplelist = ', '.join(["'{}'".format(value) for value in sample])
        samplecond = f"({samplelist})"
    if not seqid:
        seqcond = "seq_id IN (SELECT seq_id FROM sequences)"
    else:
        seqlist= [f"seq_id GLOB '{s}'" for s in seqid]
        seqcond = " OR ".join(seqlist)

    query = f"SELECT * FROM sequences WHERE {seqcond} AND sample_id in {samplecond}"
    if seqtype != "None":
        query = f"{query} AND seq_type = '{seqtype}'"

    if showquery:
        click.echo(query)
    else:
        con = sqlite3.connect(db)
        for result in query_results(con, query):
            SeqIO.write(toSeqRecord(result), output, "fasta")
        con.close()


if __name__ == "__main__":
    cli()