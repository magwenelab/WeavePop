import pandas as pd
import io
import click
from scipy import ndimage
import numpy as np
from pathlib import Path

@click.command()
@click.option('-di', '--depth_input', help='Path to BED file with depth of each window.', type=click.Path(exists=True))
@click.option('-ri', '--repeats_input', help='Path to BED file with coordinates of repetititve sequences of reference genome.', type=click.Path(exists=True))
@click.option('-sp', '--sample_name', help='Sample name as a string.', type=str)
@click.option('-wp', '--window_size', help='Size of windows in the depth BED file.', type=int)
@click.option('-dp', '--depth_threshold', help='Threshold to define copy number variation in smoothed normalzed depth.', type=click.types.FloatRange(min=0.0))
@click.option('-co', '--cnv_output', help='Path to output table of CNV calling.', type=click.Path())
def intersect_repeats(depth_input, repeats_input, sample_name, window_size, depth_threshold, cnv_output):
    print("Merge overlapping repeats and intersect with windows.")
    intersect = $(bedtools intersect -a @(depth_input) -b @(repeats_input) -wao)

    print("Reorganize intersection.")
    df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
    header = ['bed_Accession', 'bed_Start', 'bed_End','bed_Depth', 'bed_Norm_Depth', 'bed_Smooth_Depth', 'r_Accession', 'r_Start', 'r_End', 'r_Type', 'Overlap_bp'] 
    df.columns = header

    print("Calculate overlap in base pairs.")
    df = df.drop(['r_Accession', 'r_Start', 'r_End'], axis=1)
    df['r_Type_mix'] = df.groupby(['bed_Accession', 'bed_Start', 'bed_End', 'bed_Depth', 'bed_Norm_Depth', 'bed_Smooth_Depth'])['r_Type'].transform(lambda x: ','.join(x))
    df['Overlap_bp_sum'] = df.groupby(['bed_Accession', 'bed_Start', 'bed_End', 'bed_Depth', 'bed_Norm_Depth', 'bed_Smooth_Depth'])['Overlap_bp'].transform('sum')
    df = df.drop(['r_Type', 'Overlap_bp'], axis=1)
    df = df.drop_duplicates()
    df.rename(columns={'r_Type_mix': 'r_Type', 'Overlap_bp_sum': 'Overlap_bp'}, inplace=True)    
    df = df.reset_index(drop=True)

    print("Calculate fraction of window with repetitive sequences.")
    repeats_fragments = df.copy()
    repeats_fragments['Repeat_fraction'] = (repeats_fragments['Overlap_bp'] / window_size).round(2)
    repeats_fragments.columns = repeats_fragments.columns.str.replace('bed_', '')
    repeats_fragments['Sample'] = sample_name

    print("Define copy-number of regions.")
    cnv_regions = pd.DataFrame()
    for accession in repeats_fragments['Accession'].unique():
        regions_merged = repeats_fragments[repeats_fragments['Accession'] == accession].copy()
        regions_merged.loc[:, 'CNV'] = 'HAPLOID'
        regions_merged = regions_merged.reset_index(drop=True)
        for i in range(len(regions_merged)):
            if regions_merged.loc[i, 'Smooth_Depth'] > 1 + depth_threshold:
                regions_merged.loc[i, 'CNV'] = "DUPLICATION"
            elif regions_merged.loc[i, 'Smooth_Depth'] < 1 - depth_threshold:
                regions_merged.loc[i, 'CNV'] = "DELETION"
            else:
                regions_merged.loc[i, 'CNV'] = "HAPLOID"
        regions_merged.loc[:,'Region_index'] = 1
        for i in range(1, len(regions_merged)):
            if regions_merged.loc[i, 'CNV'] == regions_merged.loc[i - 1, 'CNV']:
                regions_merged.loc[i, 'Region_index'] = regions_merged.loc[i - 1, 'Region_index']
            else:
                regions_merged.loc[i, 'Region_index'] = regions_merged.loc[i - 1, 'Region_index'] + 1
        regions = regions_merged.groupby('Region_index').agg(Accession = ('Accession', 'first'),Start=('Start', 'first'), End=('End', 'last'), Depth = ('Depth', 'median'),Norm_Depth=('Norm_Depth', 'median'), Smooth_Depth=('Smooth_Depth', 'median'), CNV=('CNV', 'first'), Overlap_bp=('Overlap_bp', 'sum')).reset_index()
        
        regions['Region_Size'] = regions['End'] - regions['Start']
        regions['Repeat_fraction'] = (regions['Overlap_bp'] / regions['Region_Size']).round(2)
        regions = regions.drop(['Region_index'], axis=1)
        cnv_regions = pd.concat([cnv_regions, regions], ignore_index=True)
    print("Join regions with copy-number variants of all chromosomes.")
    cnv_regions = cnv_regions[cnv_regions['CNV'] != 'HAPLOID']
    cnv_regions = cnv_regions.round(2)
    cnv_regions['Sample'] = sample_name

    print("Convert from 0-based to 1-based coordinates")    
    cnv_regions['Start'] = cnv_regions['Start'] + 1
    # Make directory if it doesn't exist
    output_path = Path(cnv_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cnv_regions.to_csv(output_path, sep='\t', index=False, header=True)
    
if __name__ == '__main__':
    intersect_repeats()