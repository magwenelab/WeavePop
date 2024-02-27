import pandas as pd
import io
import click

@click.command()
@click.option('-s', '--structure_path', help='Path to TSV file with coordinates of putative structural variants.', type=click.Path(exists=True))
@click.option('-r', '--repeats_path', help='Path to BED file with coordinates of repetititve sequences of reference genome.', type=click.Path(exists=True))
@click.option('-o', '--output', help='Path to output file.')
@click.option('-t', '--threshold', type=float, default=0.5, help='Fraction of structural variant window size that is allowed to be covered by repetitive sequences to call it a structural variant. Default is 0.5.')

def intersect_repeats(structure_path, repeats_path, output, threshold):
    intersect = $(tail -n +2 @(structure_path) | bedtools intersect -a stdin -b @(repeats_path) -wo)

    df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
    header = ['sv_Accession', 'sv_Start', 'sv_End', 'sv_Size', 'sv_Norm_Cov', 'sv_Smooth_Cov', 'sv_Structure', 'r_Accession', 'r_Start', 'r_End', 'r_Type', 'Overlap_bp'] 
    df.columns = header
    df = df.drop(['r_Accession', 'r_Start', 'r_End'], axis=1)

    df_wide = df.pivot_table(index=['sv_Accession', 'sv_Start', 'sv_End', 'sv_Size', 'sv_Norm_Cov', 'sv_Smooth_Cov', 'sv_Structure'], columns='r_Type', values='Overlap_bp', aggfunc='sum').reset_index()
    df_wide.fillna(0, inplace=True)
    df_wide['Repeat_bp'] = df_wide.loc[:, ~df_wide.columns.str.startswith('sv_')].sum(axis=1)
    df_wide['Repeat_fraction'] = (df_wide['Repeat_bp'] / df_wide['sv_Size']).round(2)

    df_wide.columns = df_wide.columns.str.replace('sv_', '')
    df_wide['Category'] = df_wide['Repeat_fraction'].apply(lambda x: 'Structural variant' if x <= threshold else 'Window with repeats')

    df_wide.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    intersect_repeats()
