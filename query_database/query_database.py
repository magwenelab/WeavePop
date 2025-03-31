import pandas as pd
import duckdb
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO
import sys

cwd = os.getcwd()

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

def list_effect_types(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT effect_type
        FROM effects
        ORDER BY effect_type
        """
    df = con.execute(query).fetchdf()
    result = tuple(df['effect_type'])
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

def list_descriptions(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT description
        FROM gff
        ORDER BY description
        """
    df = con.execute(query).fetchdf()
    result = tuple(df['description'])
    con.close()
    return result

def get_cnv_max_length(db):
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")

    query = f"""
        SELECT MAX("region_size") AS max_length
        FROM cnvs"""
    df = con.execute(query).fetchdf()
    result = df['max_length'].values[0] 
    con.close()
    return result

def effects(db, dataset = None, sample=None, strain=None, gene_name=None, gene_id=None, impact=None, effect_type=None, lineage=None, chromosome=None, start=None, end =None):
    if gene_name and gene_id or gene_name and chromosome or gene_id and chromosome:
        raise ValueError("Only one of Gene names, Gene IDs or Location in chromosome should be provided.")
    elif sample and strain:
        raise ValueError("Only one of Sample IDs or Strains should be provided.")
    elif (sample and lineage) or (strain and lineage):
        raise ValueError("Only one of Sample IDs, Strains or Lineage should be provided.")
    elif not (gene_name or gene_id or chromosome or start or end or sample or strain or lineage):
        raise ValueError("At least one of Gene names, Gene IDs, Location, Sample IDs, Strains or Lineage should be provided.")
    else:
        pass
    if impact and effect_type:
        raise ValueError("Only one of Impacts or Effect types should be provided.")  
    
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    if not dataset:
        query_datset = f"""
            SELECT DISTINCT dataset
            FROM metadata
            """
        dataset_df = con.execute(query_datset).fetchdf()
        dataset = dataset_df['dataset'].tolist()
        dataset = tuple(dataset)
    else:
        dataset = tuple(dataset)
        
    if strain:
        strain = tuple(strain)
        query_strain = f"""
            SELECT sample
            FROM metadata
            WHERE strain IN {strain}
            """
        sample_df = con.execute(query_strain).fetchdf()
        sample = sample_df['sample'].tolist()
        sample = tuple(sample)
    elif sample:
        sample = tuple(sample)
    else:
        pass
    
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
    
    if gene_name:
        regex_pattern = '|'.join(gene_name)
        query += f"AND regexp_matches(effects.gene_name, '{regex_pattern}') "
    if gene_id:
        regex_pattern = '|'.join(gene_id)
        query += f"AND regexp_matches(effects.gene_id, '{regex_pattern}') "
    if sample:
        query += f"AND presence.sample IN {sample} "
    if lineage:
        lineage = tuple(lineage)
        query += f"AND metadata.lineage IN {lineage}"

    if chromosome :
        chromosome = tuple(chromosome)
        query += f"AND chromosomes.chromosome IN {chromosome} "
    if start or start == 0:
        query += f"AND variants.pos >= {start} "
    if end:
        query += f"AND variants.pos <= {end} "

    if impact:
        impact = tuple(impact)
        query += f"AND effects.impact IN {impact} "
    elif effect_type:
        effect_type = tuple(effect_type)
        query += f"AND effects.effect_type IN {effect_type} "
    
    print(query)
    result = con.execute(query).fetchdf()
    result = result.drop(columns=['dataset'])
    con.close()
    return result

def sequences(db, dataset=None, seq_type='DNA', sample=None, strain=None, lineage=None, gene_id=None, gene_name=None):
    if gene_name and gene_id:
        raise ValueError("Only one of Gene names or Gene IDs should be provided.")
    else:
        pass
    if sample and strain or sample and lineage or strain and lineage:
        raise ValueError("Only one of Sample IDs, Strains or Lineage should be provided.")
    else:
        pass
    if not (gene_name or gene_id or sample or strain or lineage):
        raise ValueError("At least one of Gene names, Gene IDs, Sample IDs, Strains or Lineage should be provided.")
    else:
        pass
    
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    if not dataset:
        query_datset = f"""
            SELECT DISTINCT dataset
            FROM metadata
            """
        dataset_df = con.execute(query_datset).fetchdf()
        dataset = dataset_df['dataset'].tolist()
        dataset = tuple(dataset)
    else:
        dataset = tuple(dataset)

    if gene_name:
        gene_name = tuple(gene_name)
        query_gene_id = f"""
            SELECT *
            FROM gff
            WHERE gene_name IN {gene_name}
            """
        gene_id_df = con.execute(query_gene_id).fetchdf()
        gene_id = tuple(gene_id_df['gene_id'].unique().tolist())
        print(gene_id)
    elif gene_id:
        gene_id = tuple(gene_id)
        print(gene_id)

    if strain:
        strain = tuple(strain)
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE strain IN {strain}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'].tolist())
        print(sample)
    elif lineage:
        print(lineage)
        lineage = tuple(lineage)
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE lineage IN {lineage}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'].tolist())
        print(sample)
    elif sample:
        sample = tuple(sample)
        print(sample)

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
    print(query)
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

def get_cnv(db, dataset = None, lineage=None, sample=None, strain=None, chromosome=None, cnv=None, repeat_fraction=None, start=None, end =None, min_size=None, max_size=None):
    if (sample and strain) or (sample and lineage) or (strain and lineage):
        raise ValueError("Only one of Sample IDs, Strains or Lineage should be provided.")
    
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")

    if not dataset:
        query_datset = f"""
            SELECT DISTINCT dataset
            FROM metadata
            """
        dataset_df = con.execute(query_datset).fetchdf()
        dataset = dataset_df['dataset'].tolist()
        dataset = tuple(dataset)
    else:
        dataset = tuple(dataset)
        
    query = f"""
        SELECT metadata.strain, metadata.sample, metadata.lineage, chromosomes.chromosome, cnvs.start, cnvs."end",
            cnvs.cnv, cnvs.region_size, cnvs.repeat_fraction, cnvs.feature_id,
            metadata.dataset
        FROM cnvs
        JOIN metadata ON cnvs.sample = metadata.sample
        JOIN chromosomes ON cnvs.accession = chromosomes.accession
        WHERE metadata.dataset IN {dataset}
        """
    if lineage:
        query += f"AND metadata.lineage IN {lineage}"
    if sample:
        query += f"AND metadata.sample IN {sample}"
    if strain:
        query += f"AND metadata.strain IN {strain}"
    if chromosome:
        query += f"AND chromosomes.chromosome IN {chromosome} "
    if cnv:
        query += f"AND cnvs.cnv == '{cnv}' "
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

    print(query)
    result = con.execute(query).fetchdf()
    result = result.drop(columns=['dataset'])
    con.close()
    return result

def get_metadata(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT *
        FROM metadata
        """
    df = con.execute(query).fetchdf()
    if len(df['dataset'].unique()) == 1:
        df = df.drop(columns=['dataset'])
    else:
        cols = df.columns.tolist()
        cols.insert(0, cols.pop(cols.index('dataset')))
        df = df[cols]
        
    con.close()
    return df

def genes(db, gene_name=None, gene_id=None, chromosome=None, start=None, end=None, feature_type=None, description=None, lineage=None):
    if gene_name and gene_id or gene_name and chromosome or gene_id and chromosome or gene_name and description or gene_id and description or chromosome and description or gene_name and start or gene_id and start or gene_name and end or gene_id and end:
        raise ValueError("Only one of Gene names, Gene IDs, Description or Location should be provided.")
    elif start and not end or end and not start:
        raise ValueError("Both start and end should be provided.")
    con = duckdb.connect(database=db, read_only=True)
    
    if not lineage:
        query_lineage = f"""
            SELECT DISTINCT lineage
            FROM metadata
            """
        lineage_df = con.execute(query_lineage).fetchdf()
        lineage = lineage_df['lineage'].tolist()
        lineage = tuple(lineage)
    else:
        lineage = tuple(lineage)
    
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

    if gene_name:
        query += f"AND gene_name IN {gene_name} "
    if gene_id:
        query += f"AND gene_id IN {gene_id} "
    if description:
        query += f"AND description IN {description} "
    if chromosome:
        query += f"AND chromosome IN {chromosome} "
    if start:
        query += f"AND start >= {start}"
    if end:
        query += f"""AND "end" <= {end} """
    if feature_type:
        query += f"AND primary_tag IN {feature_type}"

    print(query)
    result = con.execute(query).fetchdf()

    con.close()
    return result
