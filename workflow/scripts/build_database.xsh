import pandas as pd
import vcf
import click
import io
import os
import duckdb
import re
import sqlite3

@click.command()
@click.option('--metadata', '-m', type=click.Path(exists=True), help='Path to the sample metadata CSV file.')
@click.option('--chrom_names', '-ch', type=click.Path(), help='Path to the chromosome names file.')
@click.option('--struc_vars', '-sv', type=click.Path(), help='Path to the structural variants file.')
@click.option('--mapq_depth', '-md', type=click.Path(), help='Path to the mapq coverage file.')
@click.option('--gff', '-g', type=click.Path(), help='Path to the GFF file.')
@click.option('--effects', '-e', type=click.Path(), help='Path to the effects file of all lineages.')
@click.option('--variants', '-v', type=click.Path(), help='Path to the variants file of all lineages.')
@click.option('--presence', '-p', type=click.Path(), help='Path to the presence file of all lineages.')
@click.option('--lofs', '-l', type=click.Path(), help='Path to the lofs file of all lineages.')
@click.option('--nmds', '-n', type=click.Path(), help='Path to the nmds file of all lineages.')
@click.option('--sequences_db', '-s', type=click.Path(), help='Path to the sequences database.')
@click.option('--lineage_column', '-c', type=str, help='Name of the column in the metadata file that contains the lineage information.')
@click.option('--output', '-o', type=click.Path(), help='Output database file.')

def build_db(metadata, chrom_names, struc_vars, mapq_depth, gff, effects, variants, presence, lofs, nmds, sequences_db, output, lineage_column):
    print("Using the following arguments:")
    print(f"metadata: {metadata}")
    print(f"chrom_names: {chrom_names}")
    print(f"struc_vars: {struc_vars}")
    print(f"mapq_depth: {mapq_depth}")
    print(f"gff: {gff}")
    print(f"effects: {effects}")
    print(f"variants: {variants}")
    print(f"presence: {presence}")
    print(f"lofs: {lofs}")
    print(f"nmds: {nmds}")
    print(f"sequences_db: {sequences_db}")
    print(f"lineage_column: {lineage_column}")
    print(f"output: {output}")
    
    print("Reading and adjusting dataset files")
    df_samples = pd.read_csv(metadata)
    df_samples.columns = df_samples.columns.str.lower()
    df_samples.columns = df_samples.columns.str.replace(' ', '_')
    df_samples.rename(columns={lineage_column: 'lineage'}, inplace=True)
    print("Metadata table done!")

    df_sv = pd.read_csv(struc_vars, sep='\t')
    df_sv.columns = df_sv.columns.str.lower()
    print("Structural variants table done!")

    df_mapq_depth = pd.read_csv(mapq_depth, sep='\t')
    df_mapq_depth.rename(columns={'ID' : 'feature_id'}, inplace=True)
    print("MapQ coverage table done!")

    df_chroms = pd.read_csv(chrom_names, names=["lineage", "accession", "chromosome"], dtype = str)
    print("Chromosome names table done!")

    df_gff = pd.read_csv(gff, sep='\t', header = 0, low_memory=False)
    df_gff['feature_id_lineage'] = df_gff['feature_id'] + '_' + df_gff['lineage']
    df_gff.columns = df_gff.columns.str.lower()
    print("GFF table done!")
    gff_ids = df_gff[['feature_id','gene_id']].copy()
    gff_ids.rename(columns={'feature_id': 'transcript_id'}, inplace=True)
    gff_ids.drop_duplicates(subset='transcript_id', keep='first', inplace=True)
    print("Effects table joined with GFF IDs!")
    df_effects = pd.read_csv(effects, header = 0, sep='\t')
    df_effects = df_effects.merge(gff_ids, how='left', on='transcript_id')
    print("Effects table done!")

    df_variants = pd.read_csv(variants, header = 0, sep='\t')
    print("Variants table done!")

    df_presence = pd.read_csv(presence, header = 0, sep='\t')
    print("Presence table done!")

    df_lofs = pd.read_csv(lofs, header = 0, sep='\t')
    print("Loss of function table done!")

    df_nmds = pd.read_csv(nmds, header = 0, sep='\t')
    print("Nonsense-mediated decay table done!")

    print("Formatting dataframes done!")

    print("Connecting to database")
    con = duckdb.connect(database=output)

    print("Registering dataframes")
    con.register('df_samples', df_samples)
    con.register('df_sv', df_sv)
    con.register('df_mapq_depth', df_mapq_depth)
    con.register('df_chroms', df_chroms)
    con.register('df_gff', df_gff)
    con.register('df_effects', df_effects)
    con.register('df_variants', df_variants)
    con.register('df_lofs', df_lofs)
    con.register('df_nmds', df_nmds)
    
    print("Adding dataframes to database")
    con.execute("CREATE TABLE IF NOT EXISTS samples AS SELECT * FROM df_samples")   
    con.execute("CREATE TABLE IF NOT EXISTS structural_variants AS SELECT * FROM df_sv")
    con.execute("CREATE TABLE IF NOT EXISTS mapq_coverage AS SELECT * FROM df_mapq_depth")
    con.execute("CREATE TABLE IF NOT EXISTS chromosome_names AS SELECT * FROM df_chroms")
    con.execute("CREATE TABLE IF NOT EXISTS gff AS SELECT * FROM df_gff")
    con.execute("CREATE TABLE IF NOT EXISTS presence AS SELECT * FROM df_presence")
    con.execute("CREATE TABLE IF NOT EXISTS variants AS SELECT * FROM df_variants")
    con.execute("CREATE TABLE IF NOT EXISTS effects AS SELECT * FROM df_effects")
    con.execute("CREATE TABLE IF NOT EXISTS lofs AS SELECT * FROM df_lofs")
    con.execute("CREATE TABLE IF NOT EXISTS nmds AS SELECT * FROM df_nmds")

    print("Adding sequences database")
    seq_con = sqlite3.connect(sequences_db)
    df_seqs = pd.read_sql_query("SELECT * FROM sequences", seq_con)
    seq_con.close()
    con.register('df_seqs', df_seqs)
    con.execute("CREATE TABLE IF NOT EXISTS sequences AS SELECT * FROM df_seqs")

    print("Closing connection to database")
    con.close()
    print("Done!")

if __name__ == '__main__':
    build_db()