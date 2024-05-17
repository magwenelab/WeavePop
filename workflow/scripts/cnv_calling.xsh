import pandas as pd
import io
import click
from scipy import ndimage
import numpy as np

@click.command()
@click.option('-ci', '--coverage_input', help='Path to BED file with coverage of each region.', type=click.Path(exists=True))
@click.option('-ri', '--repeats_input', help='Path to BED file with coordinates of repetititve sequences of reference genome.', type=click.Path(exists=True))
@click.option('-chi', '--chromosome_input', help='Path to TSV file with chromosome and global modes.', type=click.Path(exists=True))
@click.option('-co', '--chromosome_output', help='Path to output table of stats by chromosome.', type=click.Path())
@click.option('-ro', '--regions_output', help='Path to output table of stats by region.', type=click.Path())
@click.option('-so', '--structural_variants_output', help='Path to output table of structural variants.', type=click.Path())
@click.option('-np', '--sample_name', help='Sample name as a string.', type=str)
@click.option('-rp', '--region_size', help='Size of regions in the coverage BED file.', type=int)
@click.option('-sp', '--smooth_size', help='Size parameter for the smoothing function.', type=int)
@click.option('-tp', '--repeats_threshold', help='Threshold for filtering out windows with too many repeats.', type=click.types.FloatRange(min=0.0))
@click.option('-cp', '--change_threshold', help='Threshold to define if region is not HAPLOID.', type=click.types.FloatRange(min=0.0))

def intersect_repeats(coverage_input, repeats_input, chromosome_input, chromosome_output, regions_output, structural_variants_output, sample_name, region_size, smooth_size, repeats_threshold,  change_threshold):
    print("Merge overlapping regions in repeats and intersect with regions.")
    intersect = $(bedtools merge -i @(repeats_input) -c 4 -o collapse | bedtools intersect -a @(coverage_input) -b stdin -wao)

    print("Reorganize intersection.")
    df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
    header = ['bed_Accession', 'bed_Start', 'bed_End', 'bed_Coverage', 'r_Accession', 'r_Start', 'r_End', 'r_Type', 'Overlap_bp'] 
    df.columns = header
    df = df.drop(['r_Accession', 'r_Start', 'r_End'], axis=1)
    df['r_Type_mix'] = df.groupby(['bed_Accession', 'bed_Start', 'bed_End', 'bed_Coverage'])['r_Type'].transform(lambda x: ','.join(x))
    df['Overlap_bp_sum'] = df.groupby(['bed_Accession', 'bed_Start', 'bed_End', 'bed_Coverage'])['Overlap_bp'].transform('sum')
    df = df.drop(['r_Type', 'Overlap_bp'], axis=1)
    df = df.drop_duplicates()
    df.rename(columns={'r_Type_mix': 'r_Type', 'Overlap_bp_sum': 'Overlap_bp'}, inplace=True)    
    df = df.reset_index(drop=True)

    print("Calculate fraction of region with repetitive sequences.")
    regions_repeats = df.copy()
    regions_repeats['Repeat_fraction'] = (regions_repeats['Overlap_bp'] / region_size).round(2)
    regions_repeats.columns = regions_repeats.columns.str.replace('bed_', '')

    print("Filter out regions with to many repeats.")
    regions_filtered = regions_repeats[regions_repeats['Repeat_fraction'] < repeats_threshold]
    regions_filtered = regions_filtered[['Accession', 'Start', 'End', 'Coverage']]
    print("Calculate genome-wide mean and median, and chromosome mean and median")
    Global_Mean = regions_filtered['Coverage'].mean().round(2)
    Global_Median = regions_filtered['Coverage'].median().round(2)
    Chrom_Mean = regions_filtered.groupby('Accession').agg(Chrom_Mean=('Coverage', 'mean')).reset_index()
    Chrom_Median = regions_filtered.groupby('Accession').agg(Chrom_Median=('Coverage', 'median')).reset_index()
    chrom_stats = pd.merge(Chrom_Mean, Chrom_Median, on='Accession')
    chrom_stats['Global_Mean'] = Global_Mean
    chrom_stats['Global_Median'] = Global_Median
    chrom_stats['Sample'] = sample_name
    chrom_stats = chrom_stats.round(2)

    chrom_mode = pd.read_csv(chromosome_input, sep='\t')
    chrom_stats = pd.merge(chrom_stats, chrom_mode, on=['Accession', 'Sample'], how="outer")
    Global_Good_Mode = chrom_stats['Global_Good_Mode'][0]
    Global_Raw_Mode = chrom_stats['Global_Raw_Mode'][0]
    chrom_stats.to_csv(chromosome_output, sep='\t', index=False, header=True)


    print("Normalize coverage.")
    regions_norm = regions_repeats[['Accession', 'Start', 'End', 'Coverage', 'Overlap_bp','Repeat_fraction']].copy()
    regions_norm.loc[:,'Norm_Mean'] = regions_norm['Coverage'] / Global_Mean
    regions_norm.loc[:,'Norm_Median'] = regions_norm['Coverage'] / Global_Median
    regions_norm.loc[:,'Norm_Mode'] = regions_norm['Coverage'] / Global_Good_Mode
    print(regions_norm.head(12))
    print("Smooth coverage.")
    cov_array = np.array(regions_norm["Norm_Mean"])
    smoothed_array = ndimage.median_filter(cov_array, size=smooth_size)
    regions_norm.loc[:,'Smooth_Mean']=pd.Series(smoothed_array)

    cov_array = np.array(regions_norm["Norm_Median"])
    smoothed_array = ndimage.median_filter(cov_array, size=smooth_size)
    regions_norm.loc[:,'Smooth_Median']=pd.Series(smoothed_array)

    cov_array = np.array(regions_norm["Norm_Mode"])
    smoothed_array = ndimage.median_filter(cov_array, size=smooth_size)
    regions_norm.loc[:,'Smooth_Mode']=pd.Series(smoothed_array)

    regions_norm = regions_norm.round(2)
    regions_norm['Sample'] = sample_name
    regions_norm.to_csv(regions_output, sep='\t', index=False, header=True)

    print("Define structure of regions.")
    structure_windows = pd.DataFrame()
    for accession in regions_norm['Accession'].unique():
        regions_windowed = regions_norm[regions_norm['Accession'] == accession].copy()
        regions_windowed.loc[:, 'Structure'] = 'HAPLOID'
        regions_windowed = regions_windowed.reset_index(drop=True)
        for i in range(len(regions_windowed)):
            if regions_windowed.loc[i, 'Smooth_Mode'] > 1 + change_threshold:
                regions_windowed.loc[i, 'Structure'] = "DUPLICATION"
            elif regions_windowed.loc[i, 'Smooth_Mode'] < 1 - change_threshold:
                regions_windowed.loc[i, 'Structure'] = "DELETION"
            else:
                regions_windowed.loc[i, 'Structure'] = "HAPLOID"

        regions_windowed.loc[:,'Window_index'] = 1
        for i in range(1, len(regions_windowed)):
            if regions_windowed.loc[i, 'Structure'] == regions_windowed.loc[i - 1, 'Structure']:
                regions_windowed.loc[i, 'Window_index'] = regions_windowed.loc[i - 1, 'Window_index']
            else:
                regions_windowed.loc[i, 'Window_index'] = regions_windowed.loc[i - 1, 'Window_index'] + 1
        windows = regions_windowed.groupby('Window_index').agg(Accession = ('Accession', 'first'),Start=('Start', 'first'), End=('End', 'last'), Coverage = ('Coverage', 'mean'),Norm_Median=('Norm_Median', 'mean'), Smooth_Mode=('Smooth_Mode', 'mean'), Structure=('Structure', 'first'), Overlap_bp=('Overlap_bp', 'sum')).reset_index()
        windows['Window_Size'] = windows['End'] - windows['Start']
        windows['Repeat_fraction'] = (windows['Overlap_bp'] / windows['Window_Size']).round(2)
        windows = windows.drop(['Window_index'], axis=1)
        structure_windows = pd.concat([structure_windows, windows], ignore_index=True)
    print("Join windows with structure variants of all chromosomes.")
    structure_windows = structure_windows[structure_windows['Structure'] != 'HAPLOID']
    structure_windows = structure_windows.round(2)
    structure_windows['Sample'] = sample_name
    structure_windows.to_csv(structural_variants_output, sep='\t', index=False, header=True)
    
if __name__ == '__main__':
    intersect_repeats()