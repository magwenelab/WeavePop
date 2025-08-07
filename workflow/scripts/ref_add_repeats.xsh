import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import pandas as pd
import io
from pathlib import Path

gff_input=snakemake.input.gff
repeats_input=snakemake.input.repeats
output=snakemake.output[0]

print("Merging annotation with repeats...")
intersect = $(bedtools intersect -a @(gff_input) -b @(repeats_input) -wao)

print("Reorganizing intersection...")
df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
header = ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attributes', 'r_Accession', 'r_Start', 'r_End', 'r_Type', 'Overlap_bp']
df.columns = header

print("Calculating overlap in base pairs...")
df = df.drop(['r_Accession', 'r_Start', 'r_End'], axis=1)
df['r_Type_mix'] = df.groupby(['attributes'])['r_Type'].transform(lambda x: ','.join(x))
df['Overlap_bp_sum'] = df.groupby(['attributes'])['Overlap_bp'].transform('sum')
df = df.drop(['r_Type', 'Overlap_bp'], axis=1)
df = df.drop_duplicates()
df.rename(columns={'r_Type_mix': 'r_Type', 'Overlap_bp_sum': 'Overlap_bp'}, inplace=True)    
df = df.reset_index(drop=True)

print("Calculating fraction of repeats...")
df['Size'] = df['end'] - df['start'] + 1
df['repeat_fraction'] = (df['Overlap_bp'] / df['Size']).round(2)
df['attributes'] = df['attributes'] + ";repeat_fraction=" + df['repeat_fraction'].astype(str)
df = df.loc[:, ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attributes']]

print("Saving df...")
filepath = Path(output)
df.to_csv(filepath, index=False, header=False, sep='\t')

print("Done!")