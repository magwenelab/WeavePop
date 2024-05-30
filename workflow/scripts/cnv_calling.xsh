import pandas as pd
import io
import click
from scipy import ndimage
import numpy as np

@click.command()
@click.option('-di', '--depth_input', help='Path to BED file with depth of each region.', type=click.Path(exists=True))
@click.option('-ri', '--repeats_input', help='Path to BED file with coordinates of repetititve sequences of reference genome.', type=click.Path(exists=True))

@click.option('-co', '--cnv_output', help='Path to output table of CNV calling.', type=click.Path())

@click.option('-sp', '--sample_name', help='Sample name as a string.', type=str)
@click.option('-rp', '--region_size', help='Size of regions in the depth BED file.', type=int)
@click.option('-dp', '--depth_threshold', help='Threshold to define copy number variation in smoothed normalzed depth.', type=click.types.FloatRange(min=0.0))

def intersect_repeats(depth_input, repeats_input, cnv_output, sample_name, region_size, depth_threshold):
    print("Merge overlapping regions in repeats and intersect with regions.")
    intersect = $(bedtools merge -i @(repeats_input) -c 4 -o collapse | bedtools intersect -a @(depth_input) -b stdin -wao)

    print("Reorganize intersection.")
    df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
    header = ['bed_Accession', 'bed_Start', 'bed_End','bed_Depth', 'bed_Norm_Depth', 'bed_Smooth_Depth', 'r_Accession', 'r_Start', 'r_End', 'r_Type', 'Overlap_bp'] 
    df.columns = header
    df = df.drop(['r_Accession', 'r_Start', 'r_End'], axis=1)
    df['r_Type_mix'] = df.groupby(['bed_Accession', 'bed_Start', 'bed_End', 'bed_Depth', 'bed_Norm_Depth', 'bed_Smooth_Depth'])['r_Type'].transform(lambda x: ','.join(x))
    df['Overlap_bp_sum'] = df.groupby(['bed_Accession', 'bed_Start', 'bed_End', 'bed_Depth', 'bed_Norm_Depth', 'bed_Smooth_Depth'])['Overlap_bp'].transform('sum')
    df = df.drop(['r_Type', 'Overlap_bp'], axis=1)
    df = df.drop_duplicates()
    df.rename(columns={'r_Type_mix': 'r_Type', 'Overlap_bp_sum': 'Overlap_bp'}, inplace=True)    
    df = df.reset_index(drop=True)

    print("Calculate fraction of region with repetitive sequences.")
    regions_repeats = df.copy()
    regions_repeats['Repeat_fraction'] = (regions_repeats['Overlap_bp'] / region_size).round(2)
    regions_repeats.columns = regions_repeats.columns.str.replace('bed_', '')
    regions_repeats['Sample'] = sample_name

    print("Define structure of regions.")
    structure_windows = pd.DataFrame()
    for accession in regions_repeats['Accession'].unique():
        regions_windowed = regions_repeats[regions_repeats['Accession'] == accession].copy()
        regions_windowed.loc[:, 'Structure'] = 'HAPLOID'
        regions_windowed = regions_windowed.reset_index(drop=True)
        for i in range(len(regions_windowed)):
            if regions_windowed.loc[i, 'Smooth_Depth'] > 1 + depth_threshold:
                regions_windowed.loc[i, 'Structure'] = "DUPLICATION"
            elif regions_windowed.loc[i, 'Smooth_Depth'] < 1 - depth_threshold:
                regions_windowed.loc[i, 'Structure'] = "DELETION"
            else:
                regions_windowed.loc[i, 'Structure'] = "HAPLOID"
        regions_windowed.loc[:,'Window_index'] = 1
        for i in range(1, len(regions_windowed)):
            if regions_windowed.loc[i, 'Structure'] == regions_windowed.loc[i - 1, 'Structure']:
                regions_windowed.loc[i, 'Window_index'] = regions_windowed.loc[i - 1, 'Window_index']
            else:
                regions_windowed.loc[i, 'Window_index'] = regions_windowed.loc[i - 1, 'Window_index'] + 1
        windows = regions_windowed.groupby('Window_index').agg(Accession = ('Accession', 'first'),Start=('Start', 'first'), End=('End', 'last'), Depth = ('Depth', 'median'),Norm_Depth=('Norm_Depth', 'median'), Smooth_Depth=('Smooth_Depth', 'median'), Structure=('Structure', 'first'), Overlap_bp=('Overlap_bp', 'sum')).reset_index()
        windows['Window_Size'] = windows['End'] - windows['Start']
        windows['Repeat_fraction'] = (windows['Overlap_bp'] / windows['Window_Size']).round(2)
        windows = windows.drop(['Window_index'], axis=1)
        structure_windows = pd.concat([structure_windows, windows], ignore_index=True)
    print("Join windows with structure variants of all chromosomes.")
    structure_windows = structure_windows[structure_windows['Structure'] != 'HAPLOID']
    structure_windows = structure_windows.round(2)
    structure_windows['Sample'] = sample_name
    structure_windows.to_csv(cnv_output, sep='\t', index=False, header=True)
    
if __name__ == '__main__':
    intersect_repeats()