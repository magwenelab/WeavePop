import pandas as pd
import duckdb
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

cwd = os.getcwd()

def list_datasets(db):
    con = duckdb.connect(database=db, read_only=True)
    
    query = f"""
            SELECT DISTINCT dataset
            FROM samples
            """
    df = con.execute(query).fetchdf()
    result = df['dataset'].tolist()
    result.sort()
    con.close()
    return result

def list_effect_types(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT effect_type
        FROM effects
        """
    df = con.execute(query).fetchdf()
    result = df['effect_type'].tolist()
    result.sort()
    result.insert(0, None)
    con.close()
    return result
    
def list_impacts(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT impact
        FROM effects
        """
    df = con.execute(query).fetchdf()
    result = df['impact'].tolist()
    result.sort()
    con.close()
    return result

def list_gene_names(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT gene_name
        FROM gff
        """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['gene_name'])
    result = df['gene_name'].tolist()
    result.sort()
    con.close()
    return result

def list_gene_ids(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT gene_id
        FROM gff
        """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['gene_id'])
    result = df['gene_id'].tolist()
    result.sort()
    con.close()
    return result

def list_samples(db, dataset=None):
    con = duckdb.connect(database=db, read_only=True)
    if dataset:
        dataset = tuple(dataset)
        query = f"""
            SELECT DISTINCT sample
            FROM samples
            WHERE dataset IN {dataset}
            """
    else:
        query = f"""
            SELECT DISTINCT sample
            FROM samples
            """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['sample'])
    result = df['sample'].tolist()
    result.sort()
    con.close()
    return result

def list_strains(db, dataset=None):
    con = duckdb.connect(database=db, read_only=True)
    if dataset:
        dataset = tuple(dataset)
        query = f"""
            SELECT DISTINCT strain
            FROM samples
            WHERE dataset IN {dataset}
            """
    else:
        query = f"""
            SELECT DISTINCT strain
            FROM samples
            """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['strain'])
    result = df['strain'].tolist()
    result.sort()
    con.close()
    return result

def list_lineages(db, dataset=None):
    con = duckdb.connect(database=db, read_only=True)
    if dataset:
        dataset = tuple(dataset)
        query = f"""
            SELECT DISTINCT lineage
            FROM samples
            WHERE dataset IN {dataset}
            """
    else:
        query = f"""
            SELECT DISTINCT lineage
            FROM samples
            """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['lineage'])
    result = df['lineage'].tolist()
    result.sort()
    result.insert(0, None)
    con.close()
    return result

def list_chromosomes(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT chromosome
        FROM chromosome_names
        """
    df = con.execute(query).fetchdf()
    result = df['chromosome'].tolist()
    result.sort()
    result.insert(0, None)
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
            FROM samples
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
            FROM samples
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
        SELECT samples.dataset, samples.strain, presence.sample, presence.lineage,
            variants.var_id, chromosome_names.chromosome,
            variants.pos AS position, variants.ref AS reference, variants.alt AS alternative,
            effects.gene_name, effects.gene_id, effects.transcript_id,
            effects.impact, effects.effect_type, effects.effect,
            effects.codon_change, effects.amino_acid_change, effects.amino_acid_length,
            effects.transcript_biotype, effects.gene_coding, effects.exon_rank,
            mapq_depth.mean_depth_normalized, mapq_depth.mean_mapq
        FROM variants 
        JOIN chromosome_names ON variants.accession = chromosome_names.accession
        JOIN presence ON variants.var_id = presence.var_id
        JOIN effects ON variants.var_id = effects.var_id
        JOIN samples ON presence.sample = samples.sample
        LEFT JOIN mapq_depth ON mapq_depth.feature_id = effects.transcript_id AND mapq_depth.sample = presence.sample
        WHERE samples.dataset IN {dataset}
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
        query += f"AND presence.lineage IN {lineage}"

    if chromosome :
        chromosome = tuple(chromosome)
        query += f"AND chromosome_names.chromosome IN {chromosome} "
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
            FROM samples
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
            FROM samples
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
            FROM samples
            WHERE lineage IN {lineage}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'].tolist())
        print(sample)
    elif sample:
        sample = tuple(sample)
        print(sample)

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

def get_cnv(db, dataset = None, lineage=None, sample=None, strain=None, chromosome=None, CNV=None, repeat_fraction=None, start=None, end =None, min_size=None, max_size=None):
    if (sample and strain) or (sample and lineage) or (strain and lineage):
        raise ValueError("Only one of Sample IDs, Strains or Lineage should be provided.")
    
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")

    if not dataset:
        query_datset = f"""
            SELECT DISTINCT dataset
            FROM samples
            """
        dataset_df = con.execute(query_datset).fetchdf()
        dataset = dataset_df['dataset'].tolist()
        dataset = tuple(dataset)
    else:
        dataset = tuple(dataset)
        
    query = f"""
        SELECT samples.strain, samples.sample, samples.lineage, chromosome_names.chromosome, cnvs.start, cnvs."end",
            cnvs.CNV, cnvs.region_size, cnvs.repeat_fraction,
            samples.dataset
        FROM cnvs
        JOIN samples ON cnvs.sample = samples.sample
        JOIN chromosome_names ON cnvs.accession = chromosome_names.accession
        WHERE samples.dataset IN {dataset}
        """
    if lineage:
        query += f"AND samples.lineage IN {lineage}"
    if sample:
        query += f"AND samples.sample IN {sample}"
    if strain:
        query += f"AND samples.strain IN {strain}"
    if chromosome:
        query += f"AND chromosome_names.chromosome IN {chromosome} "
    if CNV:
        query += f"AND cnvs.CNV == '{CNV}' "
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
        FROM samples
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

def list_feature_types(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT primary_tag
        FROM gff
        """
    df = con.execute(query).fetchdf()
    result = df['primary_tag'].tolist()
    con.close()
    return result

def list_descriptions(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT description
        FROM gff
        """
    df = con.execute(query).fetchdf()
    result = df['description'].tolist()
    result.insert(0, None)
    con.close()
    return result

def genes(db, gene_name=None, gene_id=None, chromosome=None, start=None, end=None, feature_type=None, description=None, lineage=None):
    if gene_name and gene_id or gene_name and chromosome or gene_id and chromosome or gene_name and description or gene_id and description or chromosome and description or gene_name and start or gene_id and start or gene_name and end or gene_id and end:
        raise ValueError("Only one of Gene names, Gene IDs, Description or Location should be provided.")
    elif start and not end or end and not start:
        raise ValueError("Both start and end should be provided.")
    con = duckdb.connect(database=db, read_only=True)
    
    if not lineage:
        query_lineage = f"""
            SELECT DISTINCT lineage
            FROM samples
            """
        lineage_df = con.execute(query_lineage).fetchdf()
        lineage = lineage_df['lineage'].tolist()
        lineage = tuple(lineage)
    else:
        lineage = tuple(lineage)
    
    columns = con.execute("SELECT column_name FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME = 'gff'").fetchdf()['column_name'].tolist()
    
    if 'identical_to_main_ref' and 'start_stop_mutations' in columns:
        query = f"""
            SELECT gff.lineage, chromosome_names.chromosome,
                gff.start, gff."end", gff.strand, gff.primary_tag,
                gff.gene_name, gff.gene_id,
                gff.feature_id, gff.parent,
                gff.description, gff.repeat_fraction,
                gff.identical_to_main_ref, gff.start_stop_mutations, 
            FROM gff
            JOIN chromosome_names ON gff.accession = chromosome_names.accession
            WHERE gff.lineage IN {lineage}
            """
    else:
        query = f"""
            SELECT gff.lineage, chromosome_names.chromosome,
                gff.start, gff."end", gff.strand, gff.primary_tag,
                gff.gene_name, gff.gene_id,
                gff.feature_id, gff.parent,
                gff.description, gff.repeat_fraction,
            FROM gff
            JOIN chromosome_names ON gff.accession = chromosome_names.accession
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

def mapq_depth(db, gene_name=None, gene_id=None, feature_type=None):
    if gene_name and gene_id:
        raise ValueError("Only one of Gene names or Gene IDs should be provided.")
    elif not (gene_name or gene_id):
        raise ValueError("At least one of Gene names or Gene IDs should be provided.")
    else:
        pass
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    if gene_name:
        query_gene_name = f"""
            SELECT feature_id
            FROM gff
            WHERE primary_tag = 'gene' AND gene_name IN {gene_name}
            """
        feature_id_df = con.execute(query_gene_name).fetchdf()
    elif gene_id:
        query_gene_id = f"""
            SELECT feature_id
            FROM gff
            WHERE gene_id IN {gene_id} and primary_tag = 'gene'
            """
        feature_id_df = con.execute(query_gene_id).fetchdf()
    
    feature_id = tuple(feature_id_df['feature_id'].tolist())
    
    query = f"""
        SELECT *
        FROM mapq_depth
        WHERE feature_id IN {feature_id}
        """
        
    if feature_type:
        query += f"AND primary_tag IN {feature_type}"
    print(query)
    result = con.execute(query).fetchdf()
    con.close()
    return result

