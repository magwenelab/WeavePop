import pandas as pd
import io
import click
from scipy import ndimage
import numpy as np

@click.command()
@click.option('-b', '--coverage_bam', help='Path to BED file with coverage of each region.', type=click.Path(exists=True))
@click.option('-g', '--global_input', help='Path to TSV file with chromosome and global modes.', type=click.Path(exists=True))
@click.option('-o', '--output', help='Path to output table of stats by chromosome.', type=click.Path())
@click.option('-s', '--sample_name', help='Sample name as a string.', type=str)

def intersect_repeats(coverage_input, global_input, chromosome_output, sample_name):
    
    print("Calculate chromosome mean and median")
    region_coverage = pd.read_csv(coverage_input, sep='\t')
    region_coverage.columns = ['Accession', 'Start', 'End', 'Coverage']
    Chrom_Mean = region_coverage.groupby('Accession').agg(Chrom_Mean=('Coverage', 'mean')).reset_index()
    Chrom_Median = region_coverage.groupby('Accession').agg(Chrom_Median=('Coverage', 'median')).reset_index()
    Global_Mean = region_coverage['Coverage'].mean()
    Global_Median = region_coverage['Coverage'].median()
    chrom_stats = pd.merge(Chrom_Mean, Chrom_Median, on='Accession')
    chrom_stats['Sample'] = sample_name
    chrom_stats['Global_Mean'] = Global_Mean
    chrom_stats['Global_Median'] = Global_Median
    chrom_stats = chrom_stats.round(2)
    global_mode = pd.read_csv(global_input, sep='\t')
    Global_Mode = global_mode['Global_Mode'][0]
    chrom_stats['Global_Mode'] = Global_Mode

    print("Normalize coverage")
    chrom_stats.loc[:,'Norm_Mean'] = chrom_stats['Chrom_Mean'] / chrom_stats['Global_Mode']
    chrom_stats.loc[:,'Norm_Median'] = chrom_stats['Chrom_Median'] / chrom_stats['Global_Mode']
    chrom_stats = chrom_stats.round(2)
    chrom_stats = chrom_stats[['Sample', 'Accession', 'Chrom_Mean', 'Chrom_Median', 'Norm_Mean', 'Norm_Median', 'Global_Mode', 'Global_Mean', 'Global_Median']]
    chrom_stats.to_csv(chromosome_output, sep='\t', index=False)

    
if __name__ == '__main__':
    intersect_repeats()