#/usr/bin/env xonsh -c
import pandas as pd
from pathlib import Path  
import click

@click.command()
@click.option("--mapqbed", "-mi", help='BED file with MAPQ per window', type=click.Path(exists=True))
@click.option("--depthbed", "-di", help='BED file with depth per window', type=click.Path(exists=True))
@click.option("--gff", "-gi", help='Sample gff', type=click.Path(exists=True))
@click.option("--sample", "-sp", help='Sample name', type=str)
@click.option("--depthmapq", "-dmo", help ='Intermediate BED file with MAPQ and Depth', type=click.Path())
@click.option("--output", "-o", help= 'Output table with average MAPQ and Depth per genetic feature.', type=click.Path())
def mapqdepth(mapqbed, depthbed, gff, depthmapq, sample, output):
    print("Making merged version of BED file with MAPQ and Depth columns...")
    mapq = pd.read_csv(mapqbed, names = ["chromosome", "start_w", "end_w", "mapq"], sep = "\t" )
    depth = pd.read_csv(depthbed, names = ["chromosome", "start_w", "end_w", "depth"], sep = "\t" )
    df = pd.merge(mapq, depth, on=["chromosome", "start_w", "end_w"])
    filepath = Path(depthmapq)  
    df.to_csv(filepath, index=False, header = False, sep='\t')  

    print("Intersecting MAPQ and Depth of windows with GFF file...")
    gff = $(bedtools intersect -wa -wb -a @(gff)  -b @(depthmapq))
    stringList = gff.split('\n')
    dfgff = pd.DataFrame([item.split('\t') for item in stringList], columns = ["seq_id", "source", "primary_tag", "start", "end", "score", "strand", "frame", "attribute", "chromosome", "start_w", "end_w", "mapq", "depth"])
    dfgff = dfgff.dropna()

    print("Getting average of MAPQ and Depth of all windows covered by each feature...")
    dfgff.loc[:,'mapq'] = pd.to_numeric(dfgff['mapq'], errors='coerce')
    dfgff.loc[:,'depth'] = pd.to_numeric(dfgff['depth'], errors='coerce')
    dfgff['mean_mapq'] = dfgff['mapq'].groupby(dfgff['attribute']).transform('mean')
    dfgff['mean_depth'] = dfgff['depth'].groupby(dfgff['attribute']).transform('mean')
    new_gff = dfgff.loc[:, ['seq_id', 'source', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attribute', 'mean_mapq', 'mean_depth']]
    new_gff = new_gff.drop_duplicates()
    new_gff['mean_mapq'] = new_gff['mean_mapq'].apply(lambda x: round(x, 2)).astype(str)
    new_gff['mean_depth'] = new_gff['mean_depth'].apply(lambda x: round(x, 2)).astype(str)

    print("Creating TSV with average values...")
    new_gff['feature_id'] = new_gff['attribute'].str.extract(r'ID=([^;]+)')
    tsv = new_gff.loc[:, ['feature_id','primary_tag', 'mean_mapq', 'mean_depth']]

    print("Normalizing mean_depth...")
    genome_wide_depth = df['depth'].median().round(4)
    tsv.loc[:,'mean_depth_normalized'] = tsv['mean_depth'].astype(float) / genome_wide_depth
    tsv.loc[:,'mean_depth_normalized'] = tsv.loc[:,'mean_depth_normalized'].round(2)
    
    print("Adding sample name..")
    tsv.insert(0, 'sample', sample)

    print("Saving TSV...")
    filepath = Path(output)  
    tsv.to_csv(filepath, index=False, header = True, sep='\t')

    print("Done!")
if __name__ == "__main__":
    mapqdepth()