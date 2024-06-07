#/usr/bin/env xonsh -c
import pandas as pd
from pathlib import Path  
import click

@click.command()
@click.option("--mapqbed", "-m", help='BED file with MAPQ per region', type=click.Path(exists=True))
@click.option("--covbed", "-c", help='BED file with coverage per region', type=click.Path(exists=True))
@click.option("--gff", "-g", help='Sample gff', type=click.Path(exists=True))
@click.option("--covmapq", "-cm", help ='Intermediate BED file with MAPQ and coverage', type=click.Path())
@click.option("--global_mode", "-gm", help='Global mode', type=click.Path(exists=True))
@click.option("--sample", "-s", help='Sample name', type=str)
@click.option("--output", "-o", help= 'Output table with average MAPQ and coverage per genetic feature.', type=click.Path())
def mapqcov2gff(mapqbed, covbed, gff, covmapq, global_mode, sample, output):
    print("Make merged version of BED file with MAPQ and Coverage columns")
    mapq = pd.read_csv(mapqbed, names = ["Chromosome", "Start", "End", "MAPQ"], sep = "\t" )
    cov = pd.read_csv(covbed, names = ["Chromosome", "Start", "End", "Coverage"], sep = "\t" )
    df = pd.merge(mapq, cov, on=["Chromosome", "Start", "End"])
    filepath = Path(covmapq)  
    df.to_csv(filepath, index=False, header = False, sep='\t')  

    print("Intersect MAPQ and Coverage of windows with GFF file")
    gff = $(bedtools intersect -wa -wb -a @(gff)  -b @(covmapq))
    stringList = gff.split('\n')
    dfgff = pd.DataFrame([item.split('\t') for item in stringList], columns = ["seq_id", "source", "primary_tag", "start", "end", "score", "strand", "frame", "attribute", "Chromosome", "Start", "End", "MAPQ", "COV"])
    dfgff = dfgff.dropna()

    print("Get average of MAPQ and Coverage of all windows covered by each feature")
    dfgff.loc[:,'MAPQ'] = pd.to_numeric(dfgff['MAPQ'], errors='coerce')
    dfgff.loc[:,'COV'] = pd.to_numeric(dfgff['COV'], errors='coerce')
    dfgff['mean_mapq'] = dfgff['MAPQ'].groupby(dfgff['attribute']).transform('mean')
    dfgff['mean_cov'] = dfgff['COV'].groupby(dfgff['attribute']).transform('mean')
    new_gff = dfgff.loc[:, ['seq_id', 'source', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attribute', 'mean_mapq', 'mean_cov']]
    new_gff = new_gff.drop_duplicates()
    new_gff['mean_mapq'] = new_gff['mean_mapq'].apply(lambda x: round(x, 2)).astype(str)
    new_gff['mean_cov'] = new_gff['mean_cov'].apply(lambda x: round(x, 2)).astype(str)

    print("Create TSV with average values")
    new_gff['ID'] = new_gff['attribute'].str.extract(r'ID=(.*?);')
    tsv = new_gff.loc[:, ['ID','primary_tag', 'mean_mapq', 'mean_cov']]

    print("Normalize mean_cov")
    global_df = pd.read_csv(global_mode, sep = "\t", header = 0)
    global_depth = global_df.at[0, 'Global_Mode']
    tsv['mean_cov_norm'] = tsv['mean_cov'].astype(float) / global_depth
    tsv['mean_cov_norm'] = tsv['mean_cov_norm'].round(2)
    
    print("Add sample name")
    tsv['sample'] = sample

    print("Save TSV")
    filepath = Path(output)  
    tsv.to_csv(filepath, index=False, header = True, sep='\t')

if __name__ == "__main__":
    mapqcov2gff() # Runs the function