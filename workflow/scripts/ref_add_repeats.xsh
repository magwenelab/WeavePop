import pandas as pd
import io
import click
from pathlib import Path

@click.command()
@click.option('-g', '--gff_input', help='Path to GFF file of annotation of reference.', type=click.Path(exists=True))
@click.option('-r', '--repeats_input', help='Path to BED file with coordinates of repetititve sequences of reference genome.', type=click.Path(exists=True))
@click.option('-o', '--output', help='GFF file with freaction of repeats of each genetic feature.', type=click.Path())
def intersect_repeats(gff_input, repeats_input, output):
    print("Merge annotation with repeats.")
    intersect = $(bedtools intersect -a @(gff_input) -b @(repeats_input) -wao)

    print("Reorganize intersection.")
    df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
    header = ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attributes', 'r_Accession', 'r_Start', 'r_End', 'r_Type', 'Overlap_bp']
    df.columns = header

    print("Calculate overlap in base pairs.")
    df = df.drop(['r_Accession', 'r_Start', 'r_End'], axis=1)
    df['r_Type_mix'] = df.groupby(['attributes'])['r_Type'].transform(lambda x: ','.join(x))
    df['Overlap_bp_sum'] = df.groupby(['attributes'])['Overlap_bp'].transform('sum')
    df = df.drop(['r_Type', 'Overlap_bp'], axis=1)
    df = df.drop_duplicates()
    df.rename(columns={'r_Type_mix': 'r_Type', 'Overlap_bp_sum': 'Overlap_bp'}, inplace=True)    
    df = df.reset_index(drop=True)

    print("Calculate fraction of repeats.")
    df['Size'] = df['end'] - df['start'] + 1
    df['repeat_fraction'] = (df['Overlap_bp'] / df['Size']).round(2)
    df['attributes'] = df['attributes'] + ";repeat_fraction=" + df['repeat_fraction'].astype(str)
    df = df.loc[:, ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'attributes']]

    print("Save df")
    filepath = Path(output)
    df.to_csv(filepath, index=False, header=False, sep='\t')

if __name__ == '__main__':
    intersect_repeats()