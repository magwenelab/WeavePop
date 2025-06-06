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
@click.option('--chromosomes', '-ch', type=click.Path(), help='Path to the chromosomes file.')
@click.option('--cnvs', '-cnv', type=click.Path(), help='Path to the copy-number variants file.')
@click.option('--cnv_chroms', '-cnch', type=click.Path(), help='Path to the copy-number variants per chromosome file.')
@click.option('--mapq_depth', '-md', type=click.Path(), help='Path to the mapq and depth of features file.')
@click.option('--gff_tsv', '-g', type=click.Path(), help='Path to the TSV version of the GFF file.')
@click.option('--effects', '-e', type=click.Path(), help='Path to the effects file of all lineages.')
@click.option('--variants', '-v', type=click.Path(), help='Path to the variants file of all lineages.')
@click.option('--presence', '-p', type=click.Path(), help='Path to the presence file of all lineages.')
@click.option('--lofs', '-l', type=click.Path(), help='Path to the lofs file of all lineages.')
@click.option('--nmds', '-n', type=click.Path(), help='Path to the nmds file of all lineages.')
@click.option('--sequences', '-s', type=click.Path(), help='Path to the sequences table.')
@click.option('--ref_sequences', '-r', type=click.Path(), help='Path to the reference sequences table.')
@click.option('--output', '-o', type=click.Path(), help='Output database file.')

def build_db(metadata, chromosomes, cnvs, cnv_chroms, mapq_depth, gff_tsv, effects, variants, presence, lofs, nmds, sequences, ref_sequences, output):
    print("Using the following arguments:")
    print(f"1. metadata: {metadata}")
    print(f"2. chromosomes: {chromosomes}")
    print(f"3. cnvs: {cnvs}")
    print(f"4. cnv_chroms: {cnv_chroms}")
    print(f"5. mapq_depth: {mapq_depth}")
    print(f"6. gff: {gff_tsv}")
    print(f"7. effects: {effects}")
    print(f"8. variants: {variants}")
    print(f"9. presence: {presence}")
    print(f"10. lofs: {lofs}")
    print(f"11. nmds: {nmds}")
    print(f"12. sequences: {sequences}")
    print(f"13. ref_sequences: {ref_sequences}")
    print(f"14. output: {output}")
    
    print("Reading metadata table...")
    df_metadata = pd.read_csv(metadata) 
    print("Metadata table done!")

    print("Reading Copy-number variants table...")
    df_cnv = pd.read_csv(cnvs, sep='\t')
    print("Copy-number variants table done!")

    print("Reading Copy-number variants per chromosome table...")
    df_cnv_chroms = pd.read_csv(cnv_chroms, sep='\t')
    print("Copy-number variants per chromosome table done!")

    print("Reading MAPQ-depth table...")
    df_mapq_depth = pd.read_csv(mapq_depth, sep='\t')
    print("MAPQ Depth table done!")

    print("Reading chromosome table...")
    df_chroms = pd.read_csv(chromosomes, header = 0, dtype = str)
    print("Chromosome names table done!")

    print("Reading GFF table...")
    df_gff = pd.read_csv(gff_tsv, sep='\t', header = 0, low_memory=False)
    print("GFF table done!")

    print("Reading effects table...")
    df_effects = pd.read_csv(effects, header = 0, sep='\t')
    print("Effects table done!")

    print("Reading variants table...")
    df_variants = pd.read_csv(variants, header = 0, sep='\t')
    print("Variants table done!")

    print("Reading presence table...")
    df_presence = pd.read_csv(presence, header = 0, sep='\t')
    print("Presence table done!")

    print("Reading lofs table...")
    df_lofs = pd.read_csv(lofs, header = 0, sep='\t')
    print("Loss of function table done!")

    print("Reading nonsense-mediated decay table...")
    df_nmds = pd.read_csv(nmds, header = 0, sep='\t')
    print("Nonsense-mediated decay table done!")

    print("Reading sequences table...")
    df_sequences = pd.read_csv(sequences, header = 0, sep=',')
    print("Sequences table done!")

    print("Reading reference sequences table...")
    df_ref_sequences = pd.read_csv(ref_sequences, header = 0, sep=',')
    print("Reference sequences table done!")

    print("Formatting dataframes done!")

    print("Connecting to database")
    con = duckdb.connect(database=output)

    print("Registering dataframes")
    con.register('df_metadata', df_metadata)
    con.register('df_cnv', df_cnv)
    con.register('df_cnv_chroms', df_cnv_chroms)
    con.register('df_mapq_depth', df_mapq_depth)
    con.register('df_chroms', df_chroms)
    con.register('df_gff', df_gff)
    con.register('df_effects', df_effects)
    con.register('df_variants', df_variants)
    con.register('df_lofs', df_lofs)
    con.register('df_nmds', df_nmds)
    con.register('df_sequences', df_sequences)
    con.register('df_ref_sequences', df_ref_sequences)
    
    print("Adding dataframes to database")
    con.execute("CREATE TABLE IF NOT EXISTS metadata AS SELECT * FROM df_metadata")   
    con.execute("CREATE TABLE IF NOT EXISTS cnvs AS SELECT * FROM df_cnv")
    con.execute("CREATE TABLE IF NOT EXISTS cnv_chroms AS SELECT * FROM df_cnv_chroms")
    con.execute("CREATE TABLE IF NOT EXISTS mapq_depth AS SELECT * FROM df_mapq_depth")
    con.execute("CREATE TABLE IF NOT EXISTS chromosomes AS SELECT * FROM df_chroms")
    con.execute("CREATE TABLE IF NOT EXISTS gff AS SELECT * FROM df_gff")
    con.execute("CREATE TABLE IF NOT EXISTS presence AS SELECT * FROM df_presence")
    con.execute("CREATE TABLE IF NOT EXISTS variants AS SELECT * FROM df_variants")
    con.execute("CREATE TABLE IF NOT EXISTS effects AS SELECT * FROM df_effects")
    con.execute("CREATE TABLE IF NOT EXISTS lofs AS SELECT * FROM df_lofs")
    con.execute("CREATE TABLE IF NOT EXISTS nmds AS SELECT * FROM df_nmds")
    con.execute("CREATE TABLE IF NOT EXISTS sequences AS SELECT * FROM df_sequences")
    con.execute("CREATE TABLE IF NOT EXISTS ref_sequences AS SELECT * FROM df_ref_sequences")

    print("Closing connection to database")
    con.close()
    print("Done!")

if __name__ == '__main__':
    build_db()