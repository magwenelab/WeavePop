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

def list_ref_genomes(db, dataset=None):
    con = duckdb.connect(database=db, read_only=True)
    if dataset:
        dataset = tuple(dataset)
        query = f"""
            SELECT DISTINCT ref_genome
            FROM metadata
            WHERE dataset IN {dataset}
            ORDER BY ref_genome
            """
    else:
        query = f"""
            SELECT DISTINCT ref_genome
            FROM metadata
            ORDER BY ref_genome
            """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['ref_genome'])
    result = tuple(df['ref_genome'])
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
    try:
        result_list = [int(x) for x in result_list]
    except (ValueError, TypeError):
        pass
    if all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in result_list):
        result_list.sort()
    else:
        result_list.sort(key=str)
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

def list_copy_number(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT cnv
        FROM cnv_chroms
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
        SELECT MAX("size") AS max_length
        FROM cnvs"""
    df = con.execute(query).fetchdf()
    result = df['max_length'].values[0] 
    con.close()
    return result

def effects(db, dataset = None, sample=None, strain=None, gene_name=None, gene_id=None, impact=None, effect_type=None, ref_genome=None, chromosome=None, start=None, end =None):
    if gene_name and gene_id or gene_name and chromosome or gene_id and chromosome:
        raise ValueError("Only one of Gene names, Gene IDs or Location in chromosome should be provided.")
    elif sample and strain:
        raise ValueError("Only one of Sample IDs or Strains should be provided.")
    elif (sample and ref_genome) or (strain and ref_genome):
        raise ValueError("Only one of Sample IDs, Strains or Reference genome should be provided.")
    elif not (gene_name or gene_id or chromosome or start or end or sample or strain or ref_genome):
        raise ValueError("At least one of Gene names, Gene IDs, Location, Sample IDs, Strains or Reference genome should be provided.")
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
        SELECT metadata.dataset, metadata.strain, presence.sample, metadata.ref_genome,
            variants.var_id, chromosomes.chromosome, chromosomes.accession,
            variants.pos AS position, variants.ref AS reference, variants.alt AS alternative,
            effects.gene_name, effects.gene_id, effects.feature_id,
            effects.impact, effects.effect_type, effects.effect,
            effects.codon_change, effects.amino_acid_change, effects.amino_acid_length,
            effects.transcript_biotype, effects.gene_coding, effects.exon_rank,
            variant_classification.category,
            mapq_depth.mean_depth_normalized, mapq_depth.mean_mapq
        FROM variants 
        JOIN chromosomes ON variants.accession = chromosomes.accession
        JOIN presence ON variants.var_id = presence.var_id
        JOIN effects ON variants.var_id = effects.var_id
        JOIN metadata ON presence.sample = metadata.sample
        JOIN variant_classification ON variants.var_id = variant_classification.var_id
        LEFT JOIN mapq_depth ON mapq_depth.feature_id = effects.feature_id AND mapq_depth.sample = presence.sample
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
    if ref_genome:
        ref_genome = tuple(ref_genome)
        query += f"AND metadata.ref_genome IN {ref_genome}"

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

def sequences(db, dataset=None, seq_type='DNA', sample=None, strain=None, ref_genome=None, gene_id=None, gene_name=None):
    if gene_name and gene_id:
        raise ValueError("Only one of Gene names or Gene IDs should be provided.")
    else:
        pass
    if sample and strain or sample and ref_genome or strain and ref_genome:
        raise ValueError("Only one of Sample IDs, Strains or Reference genome should be provided.")
    else:
        pass
    if not (gene_name or gene_id or sample or strain or ref_genome):
        raise ValueError("At least one of Gene names, Gene IDs, Sample IDs, Strains or Reference genome should be provided.")
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
    elif ref_genome:
        print(ref_genome)
        ref_genome = tuple(ref_genome)
        query_sample = f"""
            SELECT sample
            FROM metadata
            WHERE ref_genome IN {ref_genome}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'].tolist())
        print(sample)
    elif sample:
        sample = tuple(sample)
        print(sample)

    query = f"""
            SELECT metadata.dataset, metadata.strain, metadata.ref_genome, 
                coding_sequences.sample, coding_sequences.feature_id, coding_sequences.seq, 
                chromosomes.chromosome, chromosomes.accession,
                gff.gene_name, gff.gene_id
            FROM coding_sequences
            JOIN metadata ON coding_sequences.sample = metadata.sample
            JOIN gff ON coding_sequences.feature_id = gff.feature_id AND metadata.ref_genome = gff.ref_genome
            JOIN chromosomes ON gff.accession = chromosomes.accession
            WHERE metadata.dataset IN {dataset}"""
    if gene_id and not sample:
        query += f"""
            AND coding_sequences.feature_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )
            AND seq_type = '{seq_type}'"""
    elif gene_id and sample:
        query += f"""
            AND coding_sequences.feature_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )
            AND coding_sequences.sample IN {sample}
            AND seq_type = '{seq_type}'"""
    elif sample:
        query += f"""
            AND coding_sequences.sample IN {sample}
            AND seq_type = '{seq_type}'"""
    print(query)
    result = con.execute(query).fetchdf()
    con.close()
    return result

def ref_sequences(db, seq_type='DNA', ref_genome=None, gene_id=None, gene_name=None):
    if gene_name and gene_id:
        raise ValueError("Only one of Gene names or Gene IDs should be provided.")
    else:
        pass
    if not (gene_name or gene_id or ref_genome):
        raise ValueError("At least one of Gene names, Gene IDs, or Reference genome should be provided.")
    else:
        pass
    
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
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

    if ref_genome:
        ref_genome = tuple(ref_genome)

    query = f"""
        SELECT ref_coding_sequences.ref_genome, ref_coding_sequences.feature_id, ref_coding_sequences.seq,
                gff.gene_id, gff.gene_name,
                chromosomes.chromosome, chromosomes.accession,
        FROM ref_coding_sequences
        JOIN gff ON ref_coding_sequences.feature_id = gff.feature_id AND gff.ref_genome = ref_coding_sequences.ref_genome
        JOIN chromosomes ON gff.accession = chromosomes.accession
        WHERE seq_type = '{seq_type}'
            """
    if gene_id and not ref_genome:
        query += f"""
            AND ref_coding_sequences.feature_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )"""
    elif gene_id and ref_genome:
        query += f"""
            AND ref_coding_sequences.feature_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )
            AND ref_coding_sequences.ref_genome IN {ref_genome}
            """
    elif ref_genome:
        query += f"""
            AND ref_coding_sequences.ref_genome IN {ref_genome}
            """
    print(query)
    result = con.execute(query).fetchdf()
    con.close()
    return result

def df_to_seqrecord(df):
    records = []
    for index, row in df.iterrows():
        seq = Seq(row['seq'])
        if 'sample' in df.columns:
            record = SeqRecord(seq, id=f"{row['strain']}|{row['feature_id']}", description=f"sample={row['sample']} gene_id={row['gene_id']} gene_name={row['gene_name']} chromosome={row['chromosome']} accession={row['accession']}")
        elif 'ref_genome' in df.columns:
            record = SeqRecord(seq, id=f"{row['ref_genome']}|{row['feature_id']}", description=f"ref_genome={row['ref_genome']} gene_id={row['gene_id']} gene_name={row['gene_name']} chromosome={row['chromosome']} accession={row['accession']}")
        records.append(record)
    return records

def seqrecord_to_text(records):
    output = StringIO()
    SeqIO.write(records, output, "fasta")
    fasta_text = output.getvalue()
    output.close()
    return fasta_text

def get_cnv(db, dataset = None, ref_genome=None, sample=None, strain=None, chromosome=None, cnv=None, repeat_fraction=None, start=None, end=None, min_size=None, max_size=None):
    if (sample and strain) or (sample and ref_genome) or (strain and ref_genome):
        raise ValueError("Only one of Sample IDs, Strains or Reference genome should be provided.")
    
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
        SELECT metadata.strain, metadata.sample, metadata.ref_genome, 
            chromosomes.chromosome, chromosomes.accession,
            cnvs.start, cnvs."end",
            cnvs.size, cnvs.cnv, cnvs.depth, cnvs.norm_depth, cnvs.smooth_depth, cnvs.repeat_fraction, cnvs.repeat_overlap_bp, cnvs.feature_id,
            metadata.dataset
        FROM cnvs
        JOIN metadata ON cnvs.sample = metadata.sample
        JOIN chromosomes ON cnvs.accession = chromosomes.accession
        WHERE metadata.dataset IN {dataset}
        """
    if ref_genome:
        query += f"AND metadata.ref_genome IN {ref_genome}"
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
        query += f"AND cnvs.size >= {min_size} "
    if max_size:
        query += f"AND cnvs.size <= {max_size} "
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

def genes(db, gene_name=None, gene_id=None, chromosome=None, start=None, end=None, feature_type=None, description=None, ref_genome=None):
    if gene_name and gene_id or gene_name and chromosome or gene_id and chromosome or gene_name and description or gene_id and description or chromosome and description or gene_name and start or gene_id and start or gene_name and end or gene_id and end:
        raise ValueError("Only one of Gene names, Gene IDs, Description or Location should be provided.")
    elif start and not end or end and not start:
        raise ValueError("Both start and end should be provided.")
    con = duckdb.connect(database=db, read_only=True)
    
    if not ref_genome:
        query_ref_genome = f"""
            SELECT DISTINCT ref_genome
            FROM metadata
            """
        ref_genome_df = con.execute(query_ref_genome).fetchdf()
        ref_genome = ref_genome_df['ref_genome'].tolist()
        ref_genome = tuple(ref_genome)
    else:
        ref_genome = tuple(ref_genome)
    
    columns = con.execute("SELECT column_name FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME = 'gff'").fetchdf()['column_name'].tolist()
    
    if 'ref_identical_to_main_ref' and 'ref_start_stop_mutations' in columns:
        query = f"""
            SELECT gff.ref_genome, chromosomes.chromosome,
                chromosomes.accession,
                gff.start, gff."end", gff.strand, gff.primary_tag,
                gff.gene_name, gff.gene_id,
                gff.feature_id, gff.parent,
                gff.description, gff.ref_repeat_fraction,
                gff.ref_identical_to_main_ref, gff.ref_start_stop_mutations, 
            FROM gff
            JOIN chromosomes ON gff.accession = chromosomes.accession
            WHERE gff.ref_genome IN {ref_genome}
            """
    else:
        query = f"""
            SELECT gff.ref_genome, chromosomes.chromosome,
                chromosomes.accession,
                gff.start, gff."end", gff.strand, gff.primary_tag,
                gff.gene_name, gff.gene_id,
                gff.feature_id, gff.parent,
                gff.description, gff.ref_repeat_fraction,
            FROM gff
            JOIN chromosomes ON gff.accession = chromosomes.accession
            WHERE gff.ref_genome IN {ref_genome}
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

def get_cnv_chroms(db, dataset = None, ref_genome=None, sample=None, strain=None, chromosome=None, cnv=None, min_coverage=None, max_coverage=None):
    if (sample and strain) or (sample and ref_genome) or (strain and ref_genome):
        raise ValueError("Only one of Sample IDs, Strains or Reference genome should be provided.")
    
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
        SELECT metadata.strain, metadata.sample, metadata.ref_genome, 
            chromosomes.chromosome, chromosomes.accession, chromosomes.length,
            cnv_chroms.cnv, cnv_chroms.n_cnvs,
            cnv_chroms.total_size, cnv_chroms.coverage_percent,
            cnv_chroms.size_smallest, cnv_chroms.size_largest,
            cnv_chroms.norm_depth_mean, cnv_chroms.norm_depth_median,
            cnv_chroms.smooth_depth_mean, cnv_chroms.smooth_depth_median,
            cnv_chroms.chrom_depth, cnv_chroms.chrom_norm_depth,
             cnv_chroms.genome_depth,
            metadata.dataset
        FROM cnv_chroms
        JOIN chromosomes ON cnv_chroms.accession = chromosomes.accession
        JOIN metadata ON cnv_chroms.sample = metadata.sample
        WHERE metadata.dataset IN {dataset}
        """
    if ref_genome:
        query += f"AND metadata.ref_genome IN {ref_genome}"
    if sample:
        query += f"AND metadata.sample IN {sample}"
    if strain:
        query += f"AND metadata.strain IN {strain}"
    if chromosome:
        query += f"AND chromosomes.chromosome IN {chromosome} "
    if cnv:
        query += f"AND cnv_chroms.cnv IN {cnv} "
    if min_coverage:
        query += f"AND cnv_chroms.coverage_percent >= {min_coverage} "
    if max_coverage:
        query += f"AND cnv_chroms.coverage_percent <= {max_coverage} "
    
    print(query)
    result = con.execute(query).fetchdf()
    result = result.drop(columns=['dataset'])
    con.close()
    return result