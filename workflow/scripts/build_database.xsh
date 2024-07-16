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
@click.option('--mapq_depth', '-md', type=click.Path(), help='Path to the mapq and depth of features file.')
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
    print("MapQ Depth table done!")

    df_chroms = pd.read_csv(chrom_names, names=["lineage", "accession", "chromosome"], dtype = str)
    print("Chromosome names table done!")

    df_gff = pd.read_csv(gff, sep='\t', header = 0, low_memory=False)
    df_gff['feature_id_lineage'] = df_gff['feature_id'] + '_' + df_gff['lineage']
    df_gff.columns = df_gff.columns.str.lower()
    print("GFF table done!")
    print("Formatting effects table")
    print("Getting gene IDs from GFF file")
    gff_ids = df_gff[['gene_id', 'gene_name', 'feature_id']].copy()
    print(gff_ids.head())
    gff_ids.drop_duplicates(inplace=True)
    print("Defining function to replace gene names with gene IDs")
    def replace_with_gene_id(part):
        if part in gff_ids['gene_id'].values:
            return part
        elif part in gff_ids['gene_name'].values:
            return gff_ids.loc[gff_ids['gene_name'] == part, 'gene_id'].values[0]
        return part
    print("Reading effects table")
    effects_pre = pd.read_csv(effects, header = 0, sep='\t')
    print("Subsetting effects table")
    df_gene_transcript = effects_pre[(effects_pre['gene_name'].notnull()) & (effects_pre['transcript_id'].notnull())].copy()
    df_gene_no_transcript = effects_pre[(effects_pre['gene_name'].notnull()) & (effects_pre['transcript_id'].isnull())].copy()
    df_no_gene_no_transcript = effects_pre[(effects_pre['gene_name'].isnull()) & (effects_pre['transcript_id'].isnull())].copy()
    print("Separating fused gene names")
    df_gene_no_transcript[['part1', 'part2']] = df_gene_no_transcript['gene_name'].str.split('+', expand=True)
    print("Replacing gene names with gene IDs in part1")
    df_gene_no_transcript['gene_tag_id1'] = df_gene_no_transcript['part1'].apply(replace_with_gene_id)
    print("Replacing gene names with gene IDs in part2")
    df_gene_no_transcript.loc[df_gene_no_transcript['part2'].notnull(), 'gene_tag_id2'] = df_gene_no_transcript.loc[df_gene_no_transcript['part2'].notnull(), 'part2'].apply(replace_with_gene_id)
    print("Joining part1 with part2")
    df_gene_no_transcript['gene_id'] = df_gene_no_transcript.apply(lambda row: row['gene_tag_id1'] + '+' + row['gene_tag_id2'] if pd.notnull(row['part2']) else row['gene_tag_id1'], axis=1)
    print("Removing unnecessary columns")
    df_gene_no_transcript.drop(['part1', 'part2', 'gene_tag_id1', 'gene_tag_id2'], axis=1, inplace=True)
    print("Getting unique gene IDs from GFF")
    gff_ids_unique = gff_ids.drop_duplicates(subset='feature_id', keep='first')
    print("Creating dictionary to map feature IDs to gene IDs")
    feature_to_gene = gff_ids_unique.set_index('feature_id')['gene_id']
    print("Mapping transcript IDs to gene IDs")
    df_gene_transcript['gene_id'] = df_gene_transcript['transcript_id'].map(feature_to_gene)
    print("Concatenating dataframes")
    df_effects = pd.concat([df_gene_transcript, df_gene_no_transcript, df_no_gene_no_transcript])
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
    con.execute("CREATE TABLE IF NOT EXISTS mapq_depth AS SELECT * FROM df_mapq_depth")
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