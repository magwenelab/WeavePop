import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

from pathlib import Path
import pandas as pd

input=snakemake.input
output=snakemake.output[0]

print("Create dictionary with ref_genome and path...")
paths = {}
for path in input:
    ref_genome = Path(Path(path).stem).stem
    paths[ref_genome] = path   
    print("Available ref_genomes and paths:") 
    print(paths)
print("Files to concatenate:")    
print(paths)
print("Read files and concatenate dataframes into list...")
dfs = []
for ref_genome, path in paths.items():
    df = pd.read_csv(path, sep='\t', low_memory=False, header=0)
    df['ref_genome'] = ref_genome
    dfs.append(df)
print("Concatenate dataframes into one...")
df = pd.concat(dfs)
keep_columns = [
    'accession', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 
    'feature_id', 'parent', 'gene_id', 'gene_name', 'description', 'old_feature_id', 'ref_repeat_fraction', 'ref_genome',
    'ref_identical_to_main_ref','ref_start_stop_mutations']
existing_columns = [column for column in keep_columns if column in df.columns]
print("Keeping columns:")
print(existing_columns)
df = df[existing_columns]

print("Saving...")
df.to_csv(output, sep='\t', index=False)
print("Done!")
