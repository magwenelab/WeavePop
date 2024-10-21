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
@click.option('--cnvs', '-cnv', type=click.Path(), help='Path to the copy-number variants file.')
@click.option('--mapq_depth', '-md', type=click.Path(), help='Path to the mapq and depth of features file.')
@click.option('--gff', '-g', type=click.Path(), help='Path to the GFF file.')
@click.option('--effects', '-e', type=click.Path(), help='Path to the effects file of all lineages.')
@click.option('--variants', '-v', type=click.Path(), help='Path to the variants file of all lineages.')
@click.option('--presence', '-p', type=click.Path(), help='Path to the presence file of all lineages.')
@click.option('--lofs', '-l', type=click.Path(), help='Path to the lofs file of all lineages.')
@click.option('--nmds', '-n', type=click.Path(), help='Path to the nmds file of all lineages.')
@click.option('--sequences', '-s', type=click.Path(), help='Path to the sequences table.')
@click.option('--output', '-o', type=click.Path(), help='Output database file.')

def build_db(metadata, chrom_names, cnvs, mapq_depth, gff, effects, variants, presence, lofs, nmds, sequences, output):
    print("Using the following arguments:")
    print(f"1. metadata: {metadata}")
    print(f"2. chrom_names: {chrom_names}")
    print(f"3. cnvs: {cnvs}")
    print(f"4. mapq_depth: {mapq_depth}")
    print(f"5. gff: {gff}")
    print(f"6. effects: {effects}")
    print(f"7. variants: {variants}")
    print(f"8. presence: {presence}")
    print(f"9. lofs: {lofs}")
    print(f"10. nmds: {nmds}")
    print(f"11. sequences: {sequences}")
    print(f"12. output: {output}")
    
    print("Reading and adjusting dataset files")
    df_samples = pd.read_csv(metadata)
    df_samples.columns = df_samples.columns.str.lower()
    df_samples.columns = df_samples.columns.str.replace(' ', '_')
    if 'dataset' not in df_samples.columns:
        df_samples['dataset'] = "X"    
    print("Metadata table done!")

    print("Reading Copy-number variants table")
    df_cnv = pd.read_csv(cnvs, sep='\t')
    df_cnv.columns = df_cnv.columns.str.lower()
    print("Copy-number variants table done!")

    print("Reading MAPQ-depth table")
    df_mapq_depth = pd.read_csv(mapq_depth, sep='\t')
    df_mapq_depth.rename(columns={'ID' : 'feature_id'}, inplace=True)
    print("MAPQ Depth table done!")

    print("Reading chromosome names table")
    df_chroms = pd.read_csv(chrom_names, header = 0, dtype = str)
    print("Chromosome names table done!")

    print("Formatting GFF table")
    df_gff = pd.read_csv(gff, sep='\t', header = 0, low_memory=False)
    df_gff['feature_id_lineage'] = df_gff['feature_id'] + '_' + df_gff['lineage']
    df_gff.columns = df_gff.columns.str.lower()
    ref_mutations = ['missing_start_codon', 'missing_stop_codon', 'inframe_stop_codon']
    if all(column in df_gff.columns for column in ref_mutations):
        for column in ref_mutations:
            df_gff[column] = df_gff[column].apply(lambda x: column if x == 'Yes' else x)
        df_gff['start_stop_mutations'] = df_gff[ref_mutations].apply(lambda x: ', '.join(x.dropna()), axis=1)
        df_gff = df_gff.drop(columns=ref_mutations)
    else:
        df_gff['start_stop_mutations'] = None
    df_gff.rename(columns={'matches_ref_protein': 'identical_to_main_ref'}, inplace=True)
    print("GFF table done!")

    print("Formatting effects table")
    print("Getting gene IDs from GFF file")
    gff_ids = df_gff[['gene_id', 'gene_name', 'feature_id']].copy()
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

    print("Getting variant effects with fused gene names")
    df_fused_genes = df_gene_no_transcript[df_gene_no_transcript['gene_name'].str.contains('\\+')].copy()

    if df_fused_genes.shape[0] > 0:
        print("Fused gene names found")
        print("Separating fused gene names")
        df_fused_genes[['part1', 'part2']] = df_fused_genes['gene_name'].str.split('+', expand=True)
        print("Replacing gene names with gene IDs in part1")
        df_fused_genes['gene_tag_id1'] = df_fused_genes['part1'].apply(replace_with_gene_id)
        print("Replacing gene names with gene IDs in part2")
        df_fused_genes.loc[df_fused_genes['part2'].notnull(), 'gene_tag_id2'] = df_fused_genes.loc[df_fused_genes['part2'].notnull(), 'part2'].apply(replace_with_gene_id)
        print("Joining part1 with part2")
        df_fused_genes['gene_id'] = df_fused_genes.apply(lambda row: row['gene_tag_id1'] + '+' + row['gene_tag_id2'] if pd.notnull(row['part2']) else row['gene_tag_id1'], axis=1)
        print("Removing unnecessary columns")
        df_fused_genes.drop(['part1', 'part2', 'gene_tag_id1', 'gene_tag_id2'], axis=1, inplace=True)
        df_gene_no_transcript_fixed = df_fused_genes.copy()
    else:
        print("No fused gene names found")
        df_gene_no_transcript_fixed = df_gene_no_transcript.copy()
        df_gene_no_transcript_fixed['gene_id'] = df_gene_no_transcript_fixed['gene_name'].apply(replace_with_gene_id)

    print("Getting unique gene IDs from GFF")
    gff_ids_unique = gff_ids.drop_duplicates(subset='feature_id', keep='first')
    print("Creating dictionary to map feature IDs to gene IDs")
    feature_to_gene = gff_ids_unique.set_index('feature_id')['gene_id']
    print("Mapping transcript IDs to gene IDs")
    df_gene_transcript['gene_id'] = df_gene_transcript['transcript_id'].map(feature_to_gene)
    print("Concatenating dataframes")
    df_effects = pd.concat([df_gene_transcript, df_gene_no_transcript_fixed, df_no_gene_no_transcript])
    print("Effects table done!")

    print("Reading variants table")
    df_variants = pd.read_csv(variants, header = 0, sep='\t')
    print("Variants table done!")

    print("Reading presence table")
    df_presence = pd.read_csv(presence, header = 0, sep='\t')
    print("Presence table done!")

    print("Reading lofs table")
    df_lofs = pd.read_csv(lofs, header = 0, sep='\t')
    print("Loss of function table done!")

    print("Reading nonsense-mediated decay table")
    df_nmds = pd.read_csv(nmds, header = 0, sep='\t')
    print("Nonsense-mediated decay table done!")

    print("Reading sequences table")
    df_sequences = pd.read_csv(sequences, header = 0, sep=',')
    print("Sequences table done!")

    print("Formatting dataframes done!")

    print("Connecting to database")
    con = duckdb.connect(database=output)

    print("Registering dataframes")
    con.register('df_samples', df_samples)
    con.register('df_cnv', df_cnv)
    con.register('df_mapq_depth', df_mapq_depth)
    con.register('df_chroms', df_chroms)
    con.register('df_gff', df_gff)
    con.register('df_effects', df_effects)
    con.register('df_variants', df_variants)
    con.register('df_lofs', df_lofs)
    con.register('df_nmds', df_nmds)
    con.register('df_sequences', df_sequences)
    
    print("Adding dataframes to database")
    con.execute("CREATE TABLE IF NOT EXISTS samples AS SELECT * FROM df_samples")   
    con.execute("CREATE TABLE IF NOT EXISTS cnvs AS SELECT * FROM df_cnv")
    con.execute("CREATE TABLE IF NOT EXISTS mapq_depth AS SELECT * FROM df_mapq_depth")
    con.execute("CREATE TABLE IF NOT EXISTS chromosome_names AS SELECT * FROM df_chroms")
    con.execute("CREATE TABLE IF NOT EXISTS gff AS SELECT * FROM df_gff")
    con.execute("CREATE TABLE IF NOT EXISTS presence AS SELECT * FROM df_presence")
    con.execute("CREATE TABLE IF NOT EXISTS variants AS SELECT * FROM df_variants")
    con.execute("CREATE TABLE IF NOT EXISTS effects AS SELECT * FROM df_effects")
    con.execute("CREATE TABLE IF NOT EXISTS lofs AS SELECT * FROM df_lofs")
    con.execute("CREATE TABLE IF NOT EXISTS nmds AS SELECT * FROM df_nmds")
    con.execute("CREATE TABLE IF NOT EXISTS sequences AS SELECT * FROM df_sequences")

    print("Closing connection to database")
    con.close()
    print("Done!")

if __name__ == '__main__':
    build_db()