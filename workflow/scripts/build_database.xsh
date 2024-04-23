import pandas as pd
import vcf
import click
import io
import os
import duckdb
import re
import sqlite3

# Get the dataframes to put in the databse
def get_dataframes(lineage, db_name, temp_dir, vcf_files, db_dir, config):
    temp = os.path.join(temp_dir, "bcf_isec_" + lineage)
    if not os.path.exists(temp):
        os.makedirs(temp)

    sites_txt_file = temp + '/sites.txt'
    sites_vcf_path = temp + '/sites.vcf'
    ann_vcf_path = os.path.join(temp_dir,'snps.ann.vcf')
    bcftools_log_path = os.path.join(temp_dir, 'bcftools.log')
    snpeff_log_path = os.path.join(temp_dir, 'snpeff.log') 
    snpeff_html_path = os.path.join(temp_dir, 'snpeff.html')

    # Intersection
    print("Running bcftools isec")
    $(bcftools isec -p @(temp) @(vcf_files) 2> @(bcftools_log_path))

    # Convert sites.txt to VCF
    print("Converting sites.txt to VCF")
    sites_txt = pd.read_csv(sites_txt_file, sep='\t', header=None, names=['#CHROM', 'POS', 'REF', 'ALT', 'INFO'], dtype=str)
    sites_txt['var_id'] = 'var_' + lineage + '_' + (sites_txt.index + 1).astype(str)
    sites_vcf = sites_txt.copy()
    sites_vcf['INFO'] = 'var_id=' + sites_vcf['var_id'].astype(str) + ';' + 'MAT=' + sites_vcf['INFO'].astype(str)
    sites_vcf.drop(columns=['var_id'], inplace=True)
    sites_vcf.insert(sites_vcf.columns.get_loc('POS') + 1, 'ID', '.')
    sites_vcf.insert(sites_vcf.columns.get_loc('ALT') + 1, 'QUAL', '.')
    sites_vcf.insert(sites_vcf.columns.get_loc('QUAL') + 1, 'FILTER', '.')
    sites_vcf.to_csv(sites_vcf_path, sep='\t', index=False)

    #Get presence/absence of variants
    print("Adding sample names to presence column")
    sample_names = []
    for vcf_file in vcf_files:
        sample_names.append($(bcftools query -l @(vcf_file) 2>> @(bcftools_log_path)).rstrip())
    print(sample_names)
    print("Creating presence columns")
    presence_columns = sites_txt['INFO'].apply(lambda x: pd.Series(list(x)))
    presence_columns.columns = sample_names
    print("Creating presence matrix")
    presence_matrix = pd.concat([sites_txt, presence_columns], axis=1)
    presence_matrix = presence_matrix.drop(columns=['INFO', '#CHROM', 'POS', 'REF', 'ALT'])
    print("Converting matrix to dataframe")
    presence_melt = presence_matrix.melt(id_vars='var_id', var_name='sample', value_name='value')
    df_presence = presence_melt[presence_melt['value'] == '1'].copy()
    df_presence.drop(columns='value', inplace=True)
    # Run SnpEFF command
    print("Running SnpEff command:")
    snpeff_command = f"snpEff ann -v -classic -s {snpeff_html_path} {db_name} {sites_vcf_path} 1> {ann_vcf_path} 2> {snpeff_log_path}"
    print(snpeff_command)
    
    $(snpEff ann -v -classic -dataDir @(db_dir) -config @(config) -s @(snpeff_html_path) @(db_name) @(sites_vcf_path) 1> @(ann_vcf_path) 2> @(snpeff_log_path))

    # Get tables from SnpEFF result
    print("Getting tables from SnpEff result")
    # Create a list to store the variant data
    data_effects = []
    data_variants = []
    data_lofs = []
    data_nmds = []
    # Iterate over the records in the VCF file
    with open(ann_vcf_path, 'r') as file:
        # Open the VCF file
        vcf_reader = vcf.Reader(file)

        # Add a unique identifier to each record
        for i, record in enumerate(vcf_reader):
            var_id = record.INFO['var_id']
            var_id = var_id[0]
            data_variants.append([var_id, record.CHROM, record.POS, record.REF, record.ALT])
                    
            # Get the SnpEff annotations from the INFO field
            annotations = record.INFO.get('EFF', [])
            # Get the LOF annotations from the INFO field
            lofs = record.INFO.get('LOF', [])
            # Get the NMD annotations from the INFO field
            nmds = record.INFO.get('NMD', [])
            # Iterate over all annotations
            for annotation in annotations:
                # Split the annotation into its components
                annotation_parts = annotation.split('(')
                
                effect_type = annotation_parts[0]
                effect_info = annotation_parts[1]
                effect_info_parts = effect_info.split('|')
                impact = effect_info_parts[0]
                effect = effect_info_parts[1]
                codon_change = effect_info_parts[2]
                aa_change = effect_info_parts[3]
                aa_length = effect_info_parts[4]
                locus = effect_info_parts[5]
                biotype = effect_info_parts[6]
                coding = effect_info_parts[7]
                transcript_id = effect_info_parts[8]
                exon_rank = effect_info_parts[9]
                data_effects.append([var_id,  effect_type, impact, effect, codon_change, aa_change, aa_length, locus, biotype, coding, transcript_id, exon_rank])
            # Iterate over all lofs
            for lof in lofs:
                lof_parts = lof.split('|')
                first = lof_parts[0]
                locus = first.replace('(', '')
                nb_transcripts = lof_parts[2]
                last = lof_parts[3]
                percent_transcripts = last.replace(')', '')
                data_lofs.append([var_id, locus, nb_transcripts, percent_transcripts])
            # Iterate over all nmds
            for nmd in nmds:
                nmd_parts = nmd.split('|')
                first = nmd_parts[0]
                locus = first.replace('(', '')
                nb_transcripts = nmd_parts[2]
                last = nmd_parts[3]
                percent_transcripts = last.replace(')', '')
                data_nmds.append([var_id, locus, nb_transcripts, percent_transcripts])

    # Create dataframes from the data objects
    print("Creating dataframes")
    df_variants = pd.DataFrame(data_variants, columns=['var_id', 'accession', 'pos', 'ref', 'alt'])

    df_effects = pd.DataFrame(data_effects, columns=['var_id', 'effect_type', 'impact','effect', 'codon_change', 'amino_acid_change', 'amino_acid_length', 'locus', 'transcript_biotype', 'gene_coding', 'transcript_id', 'exon_rank'])
    df_effects['effect_id'] = 'eff_' + lineage + '_' + (df_effects.index + 1).astype(str)
    df_effects.insert(0, 'effect_id', df_effects.pop('effect_id'))

    df_lofs = pd.DataFrame(data_lofs, columns=['var_id', 'locus', 'nb_transcripts', 'percent_transcripts'])

    df_nmds = pd.DataFrame(data_nmds, columns=['var_id', 'locus', 'nb_transcripts', 'percent_transcripts'])

    dataframes = {}
    dataframes['df_presence'] = df_presence
    dataframes['df_variants'] = df_variants
    dataframes['df_effects'] = df_effects
    dataframes['df_lofs'] = df_lofs
    dataframes['df_nmds'] = df_nmds
    
    print('Adding lineage to dataframes')
    for df_name, df in dataframes.items():
        df['lineage'] = lineage

    return dataframes

# Get a tuple with the paths of the VCF files for a given lineage
def get_vcf_files(lineage, vcf_files, df_samples, lineage_column):
    lin_metadata = df_samples[df_samples[lineage_column] == lineage]
    lin_samples = lin_metadata['sample'].tolist()
    lin_vcf_files = []    
    for vcf_file in vcf_files:
        for lin_sample in lin_samples:
            if re.search(lin_sample, vcf_file):
                lin_vcf_files.append(vcf_file)
                break
    lin_vcf_files = tuple(lin_vcf_files)
    return lin_vcf_files

# # Process dataframes to include
def process_dataframes(gff,struc_vars, chrom_names, mapqcov):
    df_gff = pd.read_csv(gff, sep='\t', header = 0, low_memory=False)
    df_gff['feature_id_lineage'] = df_gff['ID'] + '_' + df_gff['lineage']
    df_gff.rename(columns={'seq_id': 'accession'}, inplace=True)
    df_gff.rename(columns={'ID' : 'feature_id'}, inplace=True)
    df_gff.columns = df_gff.columns.str.lower()

    df_sv = pd.read_csv(struc_vars, sep='\t')
    df_sv.columns = df_sv.columns.str.lower()

    df_mapqcov = pd.read_csv(mapqcov, sep='\t')
    df_mapqcov.rename(columns={'ID' : 'feature_id'}, inplace=True)

    df_chroms = pd.read_csv(chrom_names, names=["lineage", "accession", "chromosome"], dtype = str)

    dataframes = {}
    dataframes['df_sv'] = df_sv
    dataframes['df_mapqcov'] = df_mapqcov
    dataframes['df_chroms'] = df_chroms
    dataframes['df_gff'] = df_gff

    return dataframes

@click.group()
def cli():
    """
    Tool to annotate variants from sample VCF files with SnpEff and store the results in a database.
    """
    pass

@cli.command()
@click.option('--output', '-o', type=click.Path(), help='Output database file')
@click.option('--metadata', '-m', type=click.Path(exists=True), help='Path to the sample metadata CSV file')
@click.option('--lineage_column', '-c', type=str, help='Name of the column in the metadata file that contains the lineage information')
@click.option('--species_name', '-s', type=str, help='Name of the species')
@click.option('--temp_dir', '-t', type=str, help='Directory to store temporary files')
@click.option('--db_dir', '-d', type=click.Path(), help='Directory with SnpEff databse')
@click.option('--config', '-n', type=click.Path(), help='Path to the SnpEff config file')
@click.option('--struc_vars', '-v', type=click.Path(), help='Path to the structural variants file')
@click.option('--chrom_names', '-h', type=click.Path(), help='Path to the chromosome names file')
@click.option('--mapqcov', '-q', type=click.Path(), help='Path to the mapq coverage file')
@click.option('--gff', '-g', type=click.Path(), help='Path to the GFF file')
@click.option('--sequences_db', '-e', type=click.Path(), help='Path to the sequences database')
@click.argument('vcf_files', nargs=-1, type=click.Path(exists=True))
def annotate(metadata, lineage_column, species_name, temp_dir, db_dir,config, output, vcf_files, gff, struc_vars, chrom_names, mapqcov, sequences_db):
    print("Using the following parameters:")
    print(f"metadata: {metadata}")
    print(f"lineage_column: {lineage_column}")
    print(f"temp_dir: {temp_dir}")
    print(f"output: {output}")
    print(f"vcf_files: {vcf_files}")
    
    df_samples = pd.read_csv(metadata)
    df_samples.columns = df_samples.columns.str.lower()
    df_samples.columns = df_samples.columns.str.replace(' ', '_')

    df_presence = []
    df_variants = []
    df_effects = []
    df_lofs = []
    df_nmds = []
    lineages = df_samples[lineage_column].unique()
    for lineage in lineages:
        print(f"Processing lineage:")
        print(lineage)
        lin_vcf_files = get_vcf_files(lineage, vcf_files, df_samples, lineage_column)
        print("Using the following VCF files:")
        print(lin_vcf_files)
        lin_db_name = f"{species_name}_{lineage}"
        print("Using the following SnpEff database:")
        print(lin_db_name)
        lin_temp_dir = os.path.join(temp_dir, lineage)
        lin_dfs = get_dataframes(lineage, lin_db_name, lin_temp_dir, lin_vcf_files, db_dir, config)
        print("Joining dataframes")
        df_presence.append(lin_dfs['df_presence'])
        df_variants.append(lin_dfs['df_variants'])
        df_effects.append(lin_dfs['df_effects'])
        df_lofs.append(lin_dfs['df_lofs'])
        df_nmds.append(lin_dfs['df_nmds'])
        print("Finished processing lineage:")
        print(lineage)

    dataframes = {}
    dataframes['df_presence'] = pd.concat(df_presence)
    dataframes['df_variants'] = pd.concat(df_variants)
    dataframes['df_effects'] = pd.concat(df_effects)
    dataframes['df_lofs'] = pd.concat(df_lofs)
    dataframes['df_nmds'] = pd.concat(df_nmds)

    extra_dfs = process_dataframes(gff, struc_vars, chrom_names, mapqcov)

    df_samples.rename(columns={lineage_column: 'lineage'}, inplace=True)

    print("Adding dataframes to database")

    con = duckdb.connect(database=output)

    con.register('df_presence', dataframes['df_presence'])
    con.register('df_variants', dataframes['df_variants'])
    con.register('df_effects', dataframes['df_effects'])
    con.register('df_lofs', dataframes['df_lofs'])
    con.register('df_nmds', dataframes['df_nmds'])
    con.register('df_samples', df_samples)
    con.register('df_gff', extra_dfs['df_gff'])
    con.register('df_sv', extra_dfs['df_sv'])
    con.register('df_mapqcov', extra_dfs['df_mapqcov'])
    con.register('df_chroms', extra_dfs['df_chroms'])

    con.execute("CREATE TABLE IF NOT EXISTS presence AS SELECT * FROM df_presence")
    con.execute("CREATE TABLE IF NOT EXISTS variants AS SELECT * FROM df_variants")
    con.execute("CREATE TABLE IF NOT EXISTS effects AS SELECT * FROM df_effects")
    con.execute("CREATE TABLE IF NOT EXISTS lofs AS SELECT * FROM df_lofs")
    con.execute("CREATE TABLE IF NOT EXISTS nmds AS SELECT * FROM df_nmds")
    con.execute("CREATE TABLE IF NOT EXISTS samples AS SELECT * FROM df_samples")
    con.execute("CREATE TABLE IF NOT EXISTS gff AS SELECT * FROM df_gff")
    con.execute("CREATE TABLE IF NOT EXISTS structural_variants AS SELECT * FROM df_sv")
    con.execute("CREATE TABLE IF NOT EXISTS mapq_coverage AS SELECT * FROM df_mapqcov")
    con.execute("CREATE TABLE IF NOT EXISTS chromosome_names AS SELECT * FROM df_chroms")

    print("Adding sequences database")
    seq_con = sqlite3.connect(sequences_db)
    df_seqs = pd.read_sql_query("SELECT * FROM sequences", seq_con)
    seq_con.close()
    con.register('df_seqs', df_seqs)
    con.execute("CREATE TABLE IF NOT EXISTS sequences AS SELECT * FROM df_seqs")

    con.close()
    print("Done")

if __name__ == '__main__':
    cli()