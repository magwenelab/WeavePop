import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import io
from pathlib import Path
import pandas as pd

gff_input=snakemake.input.gff
output=snakemake.output[0]

find_repeats=snakemake.params.find_repeats

if find_repeats == True:
    print("Intersecting with repeats...")
    repeats_input=snakemake.input.repeats

    print("Merging annotation with repeats...")
    intersect = $(bedtools intersect -a @(gff_input) -b @(repeats_input) -wao)

    print("Reorganizing intersection...")
    df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
    header = ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attributes', 'r_accession', 'r_start', 'r_end', 'r_type', 'overlap_bp']
    df.columns = header

    print("Calculating overlap in base pairs...")
    df = df.drop(['r_accession', 'r_start', 'r_end', 'r_type'], axis=1)
    df['overlap_bp_sum'] = df.groupby(['attributes'])['overlap_bp'].transform('sum')
    df = df.drop(['overlap_bp'], axis=1)
    df = df.drop_duplicates()
    df.rename(columns={'overlap_bp_sum': 'overlap_bp'}, inplace=True)    
    df = df.reset_index(drop=True)
else:
    print("No repeats file provided, skipping intersection with repeats...")
    print("Reading gff...")
    df = pd.read_csv(gff_input, sep='\t', header=None, comment ='#')
    header = ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df.columns = header
    df['overlap_bp'] = 0


print("Calculating fraction of repeats...")
df['size'] = df['end'] - df['start'] + 1
df['repeat_fraction'] = (df['overlap_bp'] / df['size']).round(2)
df['attributes'] = df['attributes'] + ";repeat_fraction=" + df['repeat_fraction'].astype(str)
df = df.loc[:, ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attributes']]

print("Saving df...")
filepath = Path(output)
df.to_csv(filepath, index=False, header=False, sep='\t')

print("Done!")