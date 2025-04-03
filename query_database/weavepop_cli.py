import pandas as pd
import duckdb
import click
import os
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO

import query_database as qdb

cwd = os.getcwd()


# Functions to get valid values #
def list_datasets(db):
    con = duckdb.connect(database=db, read_only=True)
    
    query = f"""
            SELECT DISTINCT dataset
            FROM metadata
            ORDER BY dataset
            """
    df = con.execute(query).fetchdf()
    result = tuple(df['dataset'])
    con.close()
    return result

def list_gene_names(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT gene_name
        FROM gff
        ORDER BY gene_name
        """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['gene_name'])
    result = tuple(df['gene_name'])
    con.close()
    return result

def list_gene_ids(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT gene_id
        FROM gff
        ORDER BY gene_id
        """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['gene_id'])
    result = tuple(df['gene_id'])
    con.close()
    return result

def list_samples(db, dataset=None):
    con = duckdb.connect(database=db, read_only=True)
    if dataset:
        dataset = tuple(dataset)
        query = f"""
            SELECT DISTINCT sample
            FROM metadata
            WHERE dataset IN {dataset}
            ORDER BY sample
            """
    else:
        query = f"""
            SELECT DISTINCT sample
            FROM metadata
            ORDER BY sample
            """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['sample'])
    result = tuple(df['sample'])
    con.close()
    return result

def list_strains(db, dataset=None):
    con = duckdb.connect(database=db, read_only=True)
    if dataset:
        dataset = tuple(dataset)
        query = f"""
            SELECT DISTINCT strain
            FROM metadata
            WHERE dataset IN {dataset}
            ORDER BY strain
            """
    else:
        query = f"""
            SELECT DISTINCT strain
            FROM metadata
            ORDER BY strain
            """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['strain'])
    result = tuple(df['strain'])
    con.close()
    return result

def list_lineages(db, dataset=None):
    con = duckdb.connect(database=db, read_only=True)
    if dataset:
        dataset = tuple(dataset)
        query = f"""
            SELECT DISTINCT lineage
            FROM metadata
            WHERE dataset IN {dataset}
            ORDER BY lineage
            """
    else:
        query = f"""
            SELECT DISTINCT lineage
            FROM metadata
            ORDER BY lineage
            """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['lineage'])
    result = tuple(df['lineage'])
    con.close()
    return result

def list_chromosomes(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT chromosome
        FROM chromosomes
        """
    df = con.execute(query).fetchdf()
    result_list = df['chromosome'].tolist()
    result_list.sort(key=lambda x: (float(x) if x is not None else float('inf')))
    result = tuple(result_list)
    con.close()
    return result

def list_impacts(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT impact
        FROM effects
        ORDER BY impact
        """
    df = con.execute(query).fetchdf()
    result = tuple(df['impact'])
    con.close()
    return result

def list_feature_types(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT primary_tag
        FROM gff
        """
    df = con.execute(query).fetchdf()
    result = tuple(df['primary_tag'].unique())
    con.close()
    return result

def list_cnv_types(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT cnv
        FROM cnvs
        """
    df = con.execute(query).fetchdf()
    result_list = df['cnv'].tolist()
    result_list.insert(0, None)
    result = tuple(result_list)
    con.close()
    return result

# Function to validate input #
def validate_input(input, valid_values, valid_input):
    if all(d in valid_values for d in input):
        valid_input = input
    elif any(d not in valid_values for d in input):
        absent = tuple(d for d in input if d not in valid_values)
        print(f"{absent} not found in the database. Excluding.", file=sys.stderr)
        if len(absent) == len(input):
            print(f"No valid {valid_input} found. Exiting.", file=sys.stderr)
            sys.exit()
        else:
            valid_input = tuple(d for d in input if d in valid_values)
    return valid_input

# Functions for commands #
def get_sequences(db, dataset=None, seq_type='DNA', sample=None, strain=None, lineage=None, gene_id=None, gene_name=None):
    # Test valid combinations of input #
    if db is None:
        print("Database file must be provided. Exiting.", file=sys.stderr)
        sys.exit()
    if gene_name and gene_id:
        print("Only one of Gene names or Gene IDs should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    if sample and strain or sample and lineage or strain and lineage:
        print("Only one of Sample IDs, Strains or Lineage should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    
    # Connect to the database #
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    # Check if input values are valid and reformat input to be used in the query #
    # seq_type
    query_seq_type = f"""
        SELECT DISTINCT seq_type
        FROM sequences
        """
    seq_type_df = con.execute(query_seq_type).fetchdf()
    seq_type_tuple = tuple(seq_type_df['seq_type'])
    if seq_type not in seq_type_tuple:
        print(f"seq_type must be one of {seq_type_tuple}. Exiting.", file=sys.stderr)
        sys.exit()
    
    # dataset
    dataset_tuple = list_datasets(db)
    if dataset is None:
        dataset = dataset_tuple
    else:
        dataset_input = tuple(dataset.split(','))
        dataset = validate_input(dataset_input, dataset_tuple, 'dataset')
        
    # gene_id and gene_name
    if gene_name is not None:
        gene_name_tuple = list_gene_names(db)
        gene_name_input = tuple(gene_name.split(','))
        gene_name = validate_input(gene_name_input, gene_name_tuple, 'gene_name')
        query_gene_id = f"""
            SELECT DISTINCT gene_id
            FROM gff
            WHERE gene_name IN {gene_name}
            """
        gene_id_df = con.execute(query_gene_id).fetchdf()
        gene_id = tuple(gene_id_df['gene_id'])
    elif gene_id is not None:  
        gene_id_tuple = list_gene_ids(db) 
        gene_id_input = tuple(gene_id.split(','))
        gene_id = validate_input(gene_id_input, gene_id_tuple, 'gene_id')
    
    # sample, strain and lineage
    if strain is not None:
        strain_tuple = list_strains(db)
        strain_input = tuple(strain.split(','))
        strain = validate_input(strain_input, strain_tuple, 'strain')
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE strain IN {strain}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'])
    elif lineage is not None:
        lineage_tuple = list_lineages(db)
        lineage_input = tuple(lineage.split(','))
        lineage = validate_input(lineage_input, lineage_tuple, 'lineage')
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE lineage IN {lineage}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'])
    elif sample is not None:
        sample_tuple = list_samples(db)
        sample_input = tuple(sample.split(','))
        sample = validate_input(sample_input, sample_tuple, 'sample')
   
    # Create the query #
    query = f"""
            SELECT metadata.dataset, metadata.strain, metadata.lineage, 
                sequences.sample, sequences.transcript_id, sequences.seq, 
                chromosomes.chromosome, chromosomes.accession,
                gff.gene_name, gff.gene_id
            FROM sequences
            JOIN metadata ON sequences.sample = metadata.sample
            JOIN gff ON sequences.transcript_id = gff.feature_id AND metadata.lineage = gff.lineage
            JOIN chromosomes ON gff.accession = chromosomes.accession
            WHERE metadata.dataset IN {dataset}"""
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
    
    # Execute the query #
    result = con.execute(query).fetchdf()

    # Close the connection #
    con.close()
    return result
def get_ref_sequences(db, seq_type='DNA', lineage=None, gene_id=None, gene_name=None):
    # Test valid combinations of input #
    if db is None:
        raise ValueError("Database file must be provided.")
        sys.exit()
    if gene_name and gene_id:
        print("Only one of Gene names or Gene IDs should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    
    # Connect to the database #
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    # Check if input values are valid and reformat input to be used in the query #
    # seq_type
    query_seq_type = f"""
        SELECT DISTINCT seq_type
        FROM ref_sequences
        """
    seq_type_df = con.execute(query_seq_type).fetchdf()
    seq_type_tuple = tuple(seq_type_df['seq_type'])
    if seq_type not in seq_type_tuple:
        print(f"seq_type must be one of {seq_type_tuple}. Exiting.", file=sys.stderr)
        sys.exit()
    
    # lineage
    lineage_tuple = list_lineages(db)
    if lineage is None:
        lineage = lineage_tuple
    else:
        lineage_input = tuple(lineage.split(','))
        lineage = validate_input(lineage_input, lineage_tuple, 'lineage')
    
    # gene_id and gene_name
    if gene_name is not None:
        gene_name_tuple = list_gene_names(db)
        gene_name_input = tuple(gene_name.split(','))
        gene_name = validate_input(gene_name_input, gene_name_tuple, 'gene_name')
        query_gene_id = f"""
            SELECT DISTINCT gene_id
            FROM gff
            WHERE gene_name IN {gene_name}
            """
        gene_id_df = con.execute(query_gene_id).fetchdf()
        gene_id = tuple(gene_id_df['gene_id'])
    elif gene_id is not None:
        gene_id_tuple = list_gene_ids(db) 
        gene_id_input = tuple(gene_id.split(','))
        gene_id = validate_input(gene_id_input, gene_id_tuple, 'gene_id')
        
    # Create the query #
    query = f"""
        SELECT ref_sequences.lineage, ref_sequences.transcript_id, ref_sequences.seq,
                gff.gene_id, gff.gene_name,
                chromosomes.chromosome, chromosomes.accession,
        FROM ref_sequences
        JOIN gff ON ref_sequences.transcript_id = gff.feature_id AND gff.lineage = ref_sequences.lineage
        JOIN chromosomes ON gff.accession = chromosomes.accession
        WHERE seq_type = '{seq_type}'
            """
    if gene_id and not lineage:
        query += f"""
            AND transcript_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )"""
    elif gene_id and lineage:
        query += f"""
            AND transcript_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )
            AND ref_sequences.lineage IN {lineage}
            """
    elif lineage:
        query += f"""
            AND ref_sequences.lineage IN {lineage}
            """

    # Execute the query #
    result = con.execute(query).fetchdf()
    
    # Close the connection #
    con.close()
    return result

def df_to_seqrecord(df):
    records = []
    for index, row in df.iterrows():
        seq = Seq(row['seq'])
        if 'sample' in df.columns:
            record = SeqRecord(seq, id=f"{row['strain']}|{row['transcript_id']}", description=f"sample={row['sample']} gene_id={row['gene_id']} gene_name={row['gene_name']} chromosome={row['chromosome']} accession={row['accession']}")
        elif 'lineage' in df.columns:
            record = SeqRecord(seq, id=f"{row['lineage']}|{row['transcript_id']}", description=f"lineage={row['lineage']} gene_id={row['gene_id']} gene_name={row['gene_name']} chromosome={row['chromosome']} accession={row['accession']}")
        records.append(record)
    return records

def seqrecord_to_text(records):
    output = StringIO()
    SeqIO.write(records, output, "fasta")
    fasta_text = output.getvalue()
    output.close()
    return fasta_text

def get_variants(db, dataset=None, sample=None, strain=None, gene_name=None, gene_id=None, impact=None, 
            effect_type=None, lineage=None, chromosome=None, start=None, end=None):
    # Test valid combinations of input #
    if db is None:
        print("Database file must be provided. Exiting.", file=sys.stderr)
        sys.exit()
    if gene_name and gene_id or gene_name and chromosome or gene_id and chromosome:
        print("Only one of Gene names, Gene IDs or Location in chromosome should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    if sample and strain:
        print("Only one of Sample IDs or Strains should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    if (sample and lineage) or (strain and lineage):
        print("Only one of Sample IDs, Strains or Lineage should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    if impact and effect_type:
        print("Only one of Impacts or Effect types should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    
    # Connect to the database #
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")

    # Check if input values are valid and reformat input to be used in the query #
    # dataset
    dataset_tuple = list_datasets(db)
    
    if dataset is None:
        dataset = dataset_tuple
    else:
        dataset_input = tuple(dataset.split(','))
        dataset = validate_input(dataset_input, dataset_tuple, 'dataset')
        
    # gene_id and gene_name
    if gene_name is not None:
        gene_name_tuple = list_gene_names(db)
        gene_name_input = tuple(gene_name.split(','))
        gene_name = validate_input(gene_name_input, gene_name_tuple, 'gene_name')
        query_gene_id = f"""
            SELECT DISTINCT gene_id
            FROM gff
            WHERE gene_name IN {gene_name}
            """
        gene_id_df = con.execute(query_gene_id).fetchdf()
        gene_id = tuple(gene_id_df['gene_id'])
    elif gene_id is not None:  
        gene_id_tuple = list_gene_ids(db) 
        gene_id_input = tuple(gene_id.split(','))
        gene_id = validate_input(gene_id_input, gene_id_tuple, 'gene_id')
    
    # sample, strain and lineage
    if strain is not None:
        strain_tuple = list_strains(db)
        strain_input = tuple(strain.split(','))
        strain = validate_input(strain_input, strain_tuple, 'strain')
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE strain IN {strain}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'])
    elif lineage is not None:
        lineage_tuple = list_lineages(db)
        lineage_input = tuple(lineage.split(','))
        lineage = validate_input(lineage_input, lineage_tuple, 'lineage')
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE lineage IN {lineage}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'])
    elif sample is not None:
        sample_tuple = list_samples(db)
        sample_input = tuple(sample.split(','))
        sample = validate_input(sample_input, sample_tuple, 'sample')
                
    # impact
    if impact is not None:
        impact_tuple = list_impacts(db)
        impact_input = tuple(impact.split(','))
        impact = validate_input(impact_input, impact_tuple, 'impact')
    
    #effect_type
    if effect_type is not None:
        effect_type = tuple(effect_type.split(','))
        effect_type = tuple(e.upper() for e in effect_type)
        
    # chromosome
    if chromosome is not None:
        chromosome_tuple = list_chromosomes(db)
        chromosome_input = tuple(chromosome.split(','))
        chromosome = validate_input(chromosome_input, chromosome_tuple, 'chromosome')
    
    # Create query #
    query = f"""
        SELECT metadata.dataset, metadata.strain, presence.sample, metadata.lineage,
            variants.var_id, chromosomes.chromosome,
            variants.pos AS position, variants.ref AS reference, variants.alt AS alternative,
            effects.gene_name, effects.gene_id, effects.transcript_id,
            effects.impact, effects.effect_type, effects.effect,
            effects.codon_change, effects.amino_acid_change, effects.amino_acid_length,
            effects.transcript_biotype, effects.gene_coding, effects.exon_rank,
            mapq_depth.mean_depth_normalized, mapq_depth.mean_mapq
        FROM variants 
        JOIN chromosomes ON variants.accession = chromosomes.accession
        JOIN presence ON variants.var_id = presence.var_id
        JOIN effects ON variants.var_id = effects.var_id
        JOIN metadata ON presence.sample = metadata.sample
        LEFT JOIN mapq_depth ON mapq_depth.feature_id = effects.transcript_id AND mapq_depth.sample = presence.sample
        WHERE metadata.dataset IN {dataset}
        """
    
    if gene_id:
        regex_pattern = '|'.join(gene_id)
        query += f"AND regexp_matches(effects.gene_id, '{regex_pattern}') "
    if sample:
        query += f"AND presence.sample IN {sample} "

    if impact:
        query += f"AND effects.impact IN {impact} "
    elif effect_type:
        regex_pattern = '|'.join(effect_type)
        query += f"AND regexp_matches(effects.effect_type, '{regex_pattern}') "
    if chromosome:
        query += f"AND chromosomes.chromosome IN {chromosome} "
    if start or start == 0:
        query += f"AND variants.pos >= {start} "
    if end:
        query += f"AND variants.pos <= {end} "

    result = con.execute(query).fetchdf()
    result = result.drop(columns=['dataset'])
    con.close()
    return result

def get_cnv(db, dataset = None, lineage=None, sample=None, strain=None, chromosome=None, cnv=None, repeat_fraction=None, start=None, end =None, min_size=None, max_size=None):
    # Test valid combinations of input #
    if db is None:
        print("Database file must be provided. Exiting.", file=sys.stderr)
        sys.exit()
    if sample and strain or sample and lineage or strain and lineage:
        print("Only one of Sample IDs, Strains or Lineage should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    
    # Connect to the database #
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")

    # Check if input values are valid and reformat input to be used in the query #
    # dataset
    dataset_tuple = list_datasets(db)
    
    if dataset is None:
        dataset = dataset_tuple
    else:
        dataset_input = tuple(dataset.split(','))
        dataset = validate_input(dataset_input, dataset_tuple, 'dataset')
    
    # sample, strain and lineage
    if strain is not None:
        strain_tuple = list_strains(db)
        strain_input = tuple(strain.split(','))
        strain = validate_input(strain_input, strain_tuple, 'strain')
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE strain IN {strain}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'])
    elif lineage is not None:
        lineage_tuple = list_lineages(db)
        lineage_input = tuple(lineage.split(','))
        lineage = validate_input(lineage_input, lineage_tuple, 'lineage')
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE lineage IN {lineage}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'])
    elif sample is not None:
        sample_tuple = list_samples(db)
        sample_input = tuple(sample.split(','))
        sample = validate_input(sample_input, sample_tuple, 'sample')
    
    # cnv
    if cnv is not None:
        cnv_tuple = list_cnv_types(db)
        cnv_input = tuple(cnv.split(','))
        cnv = validate_input(cnv_input, cnv_tuple, 'cnv')
    
    # chromosome
    if chromosome is not None:
        chromosome_tuple = list_chromosomes(db)
        chromosome_input = tuple(chromosome.split(','))
        chromosome = validate_input(chromosome_input, chromosome_tuple, 'chromosome')
    
    # Create query #
    query = f"""
        SELECT metadata.strain, metadata.sample, metadata.lineage, chromosomes.chromosome, cnvs.start, cnvs."end",
            cnvs.CNV, cnvs.region_size, cnvs.repeat_fraction,
            metadata.dataset
        FROM cnvs
        JOIN metadata ON cnvs.sample = metadata.sample
        JOIN chromosomes ON cnvs.accession = chromosomes.accession
        WHERE metadata.dataset IN {dataset}
        """

    if sample:
        query += f"AND metadata.sample IN {sample}"
    if chromosome:
        query += f"AND chromosomes.chromosome IN {chromosome} "
    if cnv:
        query += f"AND cnvs.cnv IN {cnv} "
    if start or start == 0:
        query += f"AND cnvs.start >= {start} "
    if end:
        query += f"""AND cnvs."end" <= {end} """
    if min_size:
        query += f"AND cnvs.region_size >= {min_size} "
    if max_size:
        query += f"AND cnvs.region_size <= {max_size} "
    if repeat_fraction or repeat_fraction == 0:
        query += f"AND cnvs.repeat_fraction <= {repeat_fraction}"

    result = con.execute(query).fetchdf()
    result = result.drop(columns=['dataset'])
    con.close()
    return result

def get_annotation(db, lineage=None, gene_name=None, gene_id=None, chromosome=None, start=None, end=None, feature_type=None, description=None):
    # Test valid combinations of input #
    if db is None:
        print("Database file must be provided. Exiting.", file=sys.stderr)
        sys.exit()
    if (gene_name and gene_id) or (gene_name and chromosome) or (gene_id and chromosome) or (gene_name and start) or (gene_id and start) or (gene_name and end) or (gene_id and end):
        print("Only one of Gene names, Gene IDs, or Location (Chromosome, Start and End) should be provided. Exiting.", file=sys.stderr)
        sys.exit()
        
    # Connect to the database #
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    # Check if input values are valid and reformat input to be used in the query #
    # lineage
    if lineage is None:
        lineage = list_lineages(db)
    else:
        lineage_input = tuple(lineage.split(','))
        lineage_tuple = list_lineages(db)
        lineage = validate_input(lineage_input, lineage_tuple, 'lineage')
    
    # gene_id and gene_name
    if gene_name is not None:
        gene_name_tuple = list_gene_names(db)
        gene_name_input = tuple(gene_name.split(','))
        gene_name = validate_input(gene_name_input, gene_name_tuple, 'gene_name')
        query_gene_id = f"""
            SELECT DISTINCT gene_id
            FROM gff
            WHERE gene_name IN {gene_name}
            """
        gene_id_df = con.execute(query_gene_id).fetchdf()
        gene_id = tuple(gene_id_df['gene_id'])
    elif gene_id is not None:  
        gene_id_tuple = list_gene_ids(db) 
        gene_id_input = tuple(gene_id.split(','))
        gene_id = validate_input(gene_id_input, gene_id_tuple, 'gene_id')
    
    # chromosome
    if chromosome is not None:
        chromosome_tuple = list_chromosomes(db)
        chromosome_input = tuple(chromosome.split(','))
        chromosome = validate_input(chromosome_input, chromosome_tuple, 'chromosome')
    
    # start and end
    if start is not None:
        try:
            start = int(start)
        except ValueError:
            print("Start position must be an integer. Exiting.", file=sys.stderr)
            sys.exit()
    if end is not None:
        try:
            end = int(end)
        except ValueError:
            print("End position must be an integer. Exiting.", file=sys.stderr)
            sys.exit()
    
    # feature_type
    if feature_type is not None:
        feature_type_tuple = list_feature_types(db)
        feature_type_input = tuple(feature_type.split(','))
        feature_type = validate_input(feature_type_input, feature_type_tuple, 'feature_type')

    # description
    if description is not None:
        description = tuple(description.split(','))

    columns = con.execute("SELECT column_name FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME = 'gff'").fetchdf()['column_name'].tolist()
    
    if 'identical_to_main_ref' and 'start_stop_mutations' in columns:
        query = f"""
            SELECT gff.lineage, chromosomes.chromosome,
                gff.start, gff."end", gff.strand, gff.primary_tag,
                gff.gene_name, gff.gene_id,
                gff.feature_id, gff.parent,
                gff.description, gff.repeat_fraction,
                gff.identical_to_main_ref, gff.start_stop_mutations, 
            FROM gff
            JOIN chromosomes ON gff.accession = chromosomes.accession
            WHERE gff.lineage IN {lineage}
            """
    else:
        query = f"""
            SELECT gff.lineage, chromosomes.chromosome,
                gff.start, gff."end", gff.strand, gff.primary_tag,
                gff.gene_name, gff.gene_id,
                gff.feature_id, gff.parent,
                gff.description, gff.repeat_fraction,
            FROM gff
            JOIN chromosomes ON gff.accession = chromosomes.accession
            WHERE gff.lineage IN {lineage}
            """

    if gene_id:
        query += f"AND gene_id IN {gene_id} "
    if chromosome:
        query += f"AND chromosome IN {chromosome} "
    if start or start == 0:
        query += f"AND start >= {start}"
    if end:
        query += f"""AND "end" <= {end} """
    if feature_type:
        query += f"AND primary_tag IN {feature_type}"
    if description is not None:
        regex_pattern = '|'.join(description)
        query += f"AND regexp_matches(gff.description, '{regex_pattern}') "

    result = con.execute(query).fetchdf()

    con.close()
    return result

def get_metadata(db, dataset=None, lineage=None, sample=None, strain=None):
    # Test valid combinations of input #
    if db is None:
        print("Database file must be provided. Exiting.", file=sys.stderr)
        sys.exit()
    if sample and strain or sample and lineage or strain and lineage:
        print("Only one of Sample IDs, Strains or Lineage should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    
    # Connect to the database #
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    # Check if input values are valid and reformat input to be used in the query #
    # dataset
    dataset_tuple = list_datasets(db)
    
    if dataset is None:
        dataset = dataset_tuple
    else:
        dataset_input = tuple(dataset.split(','))
        dataset = validate_input(dataset_input, dataset_tuple, 'dataset')
    
    # sample, strain and lineage
    if strain is not None:
        strain_tuple = list_strains(db)
        strain_input = tuple(strain.split(','))
        strain = validate_input(strain_input, strain_tuple, 'strain')
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE strain IN {strain}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'])
    elif lineage is not None:
        lineage_tuple = list_lineages(db)
        lineage_input = tuple(lineage.split(','))
        lineage = validate_input(lineage_input, lineage_tuple, 'lineage')
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE lineage IN {lineage}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'])
    elif sample is not None:
        sample_tuple = list_samples(db)
        sample_input = tuple(sample.split(','))
        sample = validate_input(sample_input, sample_tuple, 'sample')
    
    query = f"""
        SELECT *
        FROM metadata
        WHERE dataset IN {dataset}
        """
    
    if sample:
        query += f"AND sample IN {sample} "

    result = con.execute(query).fetchdf()
    con.close()
    return result

def get_query(db, query=None, query_file=None):
    # Test valid combinations of input #
    if query is None and query_file is None:
        print("No query provided. Exiting.", file=sys.stderr)
        sys.exit()
    if query and query_file:
        print("Only one of query or query_file should be provided. Exiting.", file=sys.stderr)
        sys.exit()
    
    # Read query from file #
    if query_file:
        with open(query_file, 'r') as f:
            query = f.read()

    # Connect to the database #
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")

    # Execute the query #
    try:
        result = con.execute(query).fetchdf()
        con.close()
    except Exception as e:
        print(f"Error in query: {e}", file=sys.stderr)
        sys.exit()
    
    return result

# CLI #
@click.group()
def weavepop():
    pass

@weavepop.command()
@click.option('--db', help='Path to the database file', type=click.Path(exists=True, dir_okay=False))
@click.option('--lineage', default=None, help='Comma separated list of lineage names', type=str)
@click.option('--gene_id', default=None, help='Comma separated list of gene IDs', type=str)
@click.option('--gene_name', default=None, help='Comma separated list of gene names', type=str)
@click.option('--description', default=None, help='Comma separated list of gene descriptions', type=str)
@click.option('--chromosome', default=None, help='Comma separated list of chromosome names', type=str)
@click.option('--start', default=None, help='Minimum start position in chromosome', type=int)
@click.option('--end', default=None, help='Maximum end position in chromosome', type=int)
@click.option('--feature_type', default=None, help='Comma separated list of feature types.', type=str)
@click.option('--output', default=None, help='Path to output file. Printed to standard output if not provided.', type=click.File('w'))
def annotation(db, lineage, gene_id, gene_name, description, chromosome, start, end, feature_type, output):
    result = get_annotation(db, lineage=lineage, gene_id=gene_id, gene_name=gene_name, description=description,
                                chromosome=chromosome, start=start, end=end, feature_type=feature_type)
    if output is None:
        result.to_csv(sys.stdout, sep='\t',index=False)
    else:
        result.to_csv(output, sep='\t', index=False)
    

@weavepop.command()
@click.option('--db', help='Path to the database file', type=click.Path(exists=True, dir_okay=False))
@click.option('--dataset', default=None, help='Comma separated list of dataset names', type=str)
@click.option('--gene_id', default=None, help='Comma separated list of gene IDs', type=str)
@click.option('--gene_name', default=None, help='Comma separated list of gene names', type=str)
@click.option('--sample', default=None, help='Comma separated list of sample IDs', type=str)
@click.option('--strain', default=None, help='Comma separated list of strain names', type=str)
@click.option('--lineage', default=None, help='Comma separated list of lineage names', type=str)
@click.option('--seq_type', default='DNA', help='Sequence type. Options are DNA and PROTEIN', show_default=True, type=str)
@click.option('--output', default=None, help='Path to output file. Printed to standard output if not provided.', type=click.File('w'))
def sequences(db, gene_id, gene_name, dataset, sample, strain, lineage, seq_type, output):
    result = get_sequences(db, gene_id=gene_id, gene_name=gene_name, dataset=dataset, sample=sample, 
                           strain=strain, lineage=lineage, seq_type=seq_type)
    records = df_to_seqrecord(result)
    fasta_text = seqrecord_to_text(records)
    if output is None:
        print(fasta_text)
    else:
        with open(output, 'w') as f:
            f.write(fasta_text)

@weavepop.command()
@click.option('--db', help='Path to the database file', type=click.Path(exists=True, dir_okay=False))
@click.option('--gene_id', default=None, help='Comma separated list of gene IDs', type=str)
@click.option('--gene_name', default=None, help='Comma separated list of gene names', type=str)
@click.option('--lineage', default=None, help='Comma separated list of lineage names', type=str)
@click.option('--seq_type', default='DNA', help='Sequence type. Options are DNA and PROTEIN', show_default=True, type=str)
@click.option('--output', default=None, help='Path to output file. Printed to standard output if not provided.', type=click.File('w'))
def ref_sequences(db, gene_id, gene_name,  lineage, seq_type, output):
    result = get_ref_sequences(db, gene_id=gene_id, gene_name=gene_name, 
                           lineage=lineage, seq_type=seq_type)
    records = df_to_seqrecord(result)
    fasta_text = seqrecord_to_text(records)
    if output is None:
        print(fasta_text)
    else:
        with open(output, 'w') as f:
            f.write(fasta_text)

@weavepop.command()
@click.option('--db', help='Path to the database file', type=click.Path(exists=True, dir_okay=False))
@click.option('--dataset', default=None, help='Comma separated list of dataset names', type=str)
@click.option('--gene_id', default=None, help='Comma separated list of gene IDs', type=str)
@click.option('--gene_name', default=None, help='Comma separated list of gene names', type=str)
@click.option('--sample', default=None, help='Comma separated list of sample IDs', type=str)
@click.option('--strain', default=None, help='Comma separated list of strain names', type=str)
@click.option('--lineage', default=None, help='Comma separated list of lineage names', type=str)
@click.option('--impact', default=None, help='Comma separated list of impacts', type=str)
@click.option('--effect_type', default=None, help='Comma separated list of effect types', type=str)
@click.option('--chromosome', default=None, help='Comma separated list of chromosome names', type=str)
@click.option('--start', default=None, help='Minimum start position in chromosome', type=int)
@click.option('--end', default=None, help='Maximum end position in chromosome', type=int)
@click.option('--output', default=None, help='Path to output file. Printed to standard output if not provided.', type=click.File('w'))
def variants(db, dataset, gene_id, gene_name, sample, strain, lineage, impact, effect_type, chromosome, 
             start, end, output):
    result = get_variants(db, dataset=dataset, gene_id=gene_id, gene_name=gene_name,  sample=sample, strain=strain, 
                     lineage=lineage, impact=impact, effect_type=effect_type, chromosome=chromosome, 
                     start=start, end=end)
    if output is None:
        result.to_csv(sys.stdout, sep='\t',index=False)
    else:
        result.to_csv(output, sep='\t', index=False)
        
@weavepop.command()
@click.option('--db', help='Path to the database file', type=click.Path(exists=True, dir_okay=False))
@click.option('--dataset', default=None, help='Comma separated list of dataset names', type=str)
@click.option('--strain', default=None, help='Comma separated list of strain names', type=str)
@click.option('--sample', default=None, help='Comma separated list of sample IDs', type=str)
@click.option('--lineage', default=None, help='Comma separated list of lineage names', type=str)
@click.option('--chromosome', default=None, help='Comma separated list of chromosome names', type=str)
@click.option('--start', default=None, help='Minimum start position in chromosome', type=int)
@click.option('--end', default=None, help='Maximum end position in chromosome', type=int)
@click.option('--min_size', default=None, help='Minimum size of CNV', type=int)
@click.option('--max_size', default=None, help='Maximum size of CNV', type=int)
@click.option('--cnv', default=None, help='Comma separated list of CNV types. Options are duplication and deletion.', type=str)
@click.option('--repeat_fraction', default=None, help='Max repeat fraction allowed in CNV', type=float)
@click.option('--output', default=None, help='Path to output file. Printed to standard output if not provided.', type=click.File('w'))
def cnv(db, dataset, strain, sample, lineage, chromosome, start, end, min_size, max_size, cnv, repeat_fraction, output):
    result = get_cnv(db, dataset=dataset ,strain=strain, sample=sample, lineage=lineage,
                    chromosome=chromosome, start=start, end =end, min_size=min_size, max_size=max_size,
                    cnv=cnv, repeat_fraction=repeat_fraction)
    if output is None:
        result.to_csv(sys.stdout, sep='\t',index=False)
    else:
        result.to_csv(output, sep='\t', index=False)

@weavepop.command()
@click.option('--db', help='Path to the database file', type=click.Path(exists=True, dir_okay=False))
@click.option('--dataset', default=None, help='Comma separated list of dataset names', type=str)
@click.option('--lineage', default=None, help='Comma separated list of lineage names', type=str)
@click.option('--sample', default=None, help='Comma separated list of sample IDs', type=str)
@click.option('--strain', default=None, help='Comma separated list of strain names', type=str)
@click.option('--output', default=None, help='Path to output file. Printed to standard output if not provided.', type=click.File('w'))
def metadata(db, dataset, lineage, sample, strain, output):
    result = get_metadata(db, dataset=dataset, lineage=lineage, sample=sample, strain=strain)
    if output is None:
        result.to_csv(sys.stdout, sep='\t',index=False)
    else:
        result.to_csv(output, sep='\t', index=False)

@weavepop.command()
@click.option('--db', help='Path to the database file', type=click.Path(exists=True, dir_okay=False))
@click.option('--output', default=None, help='Path to output file. Printed to standard output if not provided.', type=click.File('w'))
@click.option('--query', default=None, help='SQL query to execute', type=str)
@click.option("--query_file",default=None, help='File with SQL query to execute', type=click.Path(exists=True, dir_okay=False))
def query(db, output, query, query_file):
    result = get_query(db, query=query, query_file=query_file)
    if output is None:
        result.to_csv(sys.stdout, sep='\t',index=False)
    else:
        result.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    weavepop()
