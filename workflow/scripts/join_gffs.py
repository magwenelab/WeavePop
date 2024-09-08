import click
import pandas as pd
from pathlib import Path

@click.command()
@click.option('--output', '-o', required=True, type=click.Path(), help="Path to the output file")
@click.argument('lineage_tsv', nargs=-1, type=click.Path(exists=True))   

def join_gff(output, lineage_tsv):
    print("Create dictionary with lineage and path")
    paths = {}
    for path in lineage_tsv:
        lineage = Path(Path(path).stem).stem
        paths[lineage] = path   
        print("Available lineages and paths:") 
        print(lineage, path)
    print("Files to concatenate:")    
    print(paths)
    print("Read files and concatenate dataframes into list")
    dfs = []
    for lineage, path in paths.items():
        df = pd.read_csv(path, sep='\t', low_memory=False, header=0)
        df['lineage'] = lineage
        dfs.append(df)
    print("Concatenate dataframes into one")
    df = pd.concat(dfs)
    keep_columns = [
        'seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 
        'frame', 'ID', 'Parent', 'locus', 'Name', 'description', 'old_ID', 'lineage',
        'matches_ref_protein','missing_start_codon', 'missing_stop_codon', 'inframe_stop_codon']
    existing_columns = [column for column in keep_columns if column in df.columns]
    print(existing_columns)
    df = df[existing_columns]
    print("Replacing True/False for Yes/No")
    df = df.replace({True: 'Yes', False: 'No'})
    print("Renaming columns")
    df = df.rename(columns={
        'seq_id': 'accession',
        'ID': 'feature_id',
        'Name': 'gene_name',
        'locus': 'gene_id',
        'old_ID': 'old_feature_id'})
    print("Saving")
    df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    join_gff()   