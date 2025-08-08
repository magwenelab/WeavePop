import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

from pathlib import Path  
import pandas as pd

mapqbed=snakemake.input.mapqbed
depthbed=snakemake.input.depthbed
gff=snakemake.input.gff
depthmapq=snakemake.output.depthmapq
sample=snakemake.wildcards.sample
output=snakemake.output.tsv

print("Making merged version of BED file with MAPQ and Depth columns...")
mapq = pd.read_csv(mapqbed, names = ["chromosome", "start_w", "end_w", "mapq"], sep = "\t" )
depth = pd.read_csv(depthbed, names = ["chromosome", "start_w", "end_w", "depth"], sep = "\t" )
df = pd.merge(mapq, depth, on=["chromosome", "start_w", "end_w"])
filepath = Path(depthmapq)  
df.to_csv(filepath, index=False, header = False, sep='\t')  

print("Intersecting MAPQ and Depth of windows with GFF file...")
gff = $(bedtools intersect -wa -wb -a @(gff)  -b @(depthmapq) 2>> @(log_file))
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