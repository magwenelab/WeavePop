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
def cnv_calling(depth_input, repeats_input, sample_name, window_size, depth_threshold, cnv_output):
    print("Merge overlapping repeats and intersect with windows.")
    intersect = $(bedtools intersect -a @(depth_input) -b @(repeats_input) -wao)

    print("Reorganize intersection.")
    df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
    header = ['bed_accession', 'bed_start', 'bed_end','bed_depth', 'bed_norm_depth', 'bed_smooth_depth', 'r_accession', 'r_start', 'r_end', 'r_type', 'overlap_bp'] 
    df.columns = header

    print("Calculate overlap in base pairs.")
    df = df.drop(['r_accession', 'r_start', 'r_end'], axis=1)
    df['r_type_mix'] = df.groupby(['bed_accession', 'bed_start', 'bed_end', 'bed_depth', 'bed_norm_depth', 'bed_smooth_depth'])['r_type'].transform(lambda x: ','.join(x))
    df['overlap_bp_sum'] = df.groupby(['bed_accession', 'bed_start', 'bed_end', 'bed_depth', 'bed_norm_depth', 'bed_smooth_depth'])['overlap_bp'].transform('sum')
    df = df.drop(['r_type', 'overlap_bp'], axis=1)
    df = df.drop_duplicates()
    df.rename(columns={'r_type_mix': 'r_type', 'overlap_bp_sum': 'overlap_bp'}, inplace=True)    
    df = df.reset_index(drop=True)

    print("Calculate fraction of window with repetitive sequences.")
    repeats_fragments = df.copy()
    repeats_fragments['repeat_fraction'] = (repeats_fragments['overlap_bp'] / window_size).round(2)
    repeats_fragments.columns = repeats_fragments.columns.str.replace('bed_', '')
    repeats_fragments['sample'] = sample_name

    print("Define copy-number of regions.")
    cnv_regions = pd.DataFrame()
    for accession in repeats_fragments['accession'].unique():
        regions_merged = repeats_fragments[repeats_fragments['accession'] == accession].copy()
        regions_merged.loc[:, 'cnv'] = 'haploid'
        regions_merged = regions_merged.reset_index(drop=True)
        for i in range(len(regions_merged)):
            if regions_merged.loc[i, 'smooth_depth'] > 1 + depth_threshold:
                regions_merged.loc[i, 'cnv'] = "duplication"
            elif regions_merged.loc[i, 'smooth_depth'] < 1 - depth_threshold:
                regions_merged.loc[i, 'cnv'] = "deletion"
            else:
                regions_merged.loc[i, 'cnv'] = "haploid"
        regions_merged.loc[:,'region_index'] = 1
        for i in range(1, len(regions_merged)):
            if regions_merged.loc[i, 'cnv'] == regions_merged.loc[i - 1, 'cnv']:
                regions_merged.loc[i, 'region_index'] = regions_merged.loc[i - 1, 'region_index']
            else:
                regions_merged.loc[i, 'region_index'] = regions_merged.loc[i - 1, 'region_index'] + 1
        regions = regions_merged.groupby('region_index').agg(accession = ('accession', 'first'),start=('start', 'first'), end=('end', 'last'), depth = ('depth', 'median'),norm_depth=('norm_depth', 'median'), smooth_depth=('smooth_depth', 'median'), cnv=('cnv', 'first'), overlap_bp=('overlap_bp', 'sum')).reset_index()
        
        regions['region_size'] = regions['end'] - regions['start']
        regions['repeat_fraction'] = (regions['overlap_bp'] / regions['region_size']).round(2)
        regions = regions.drop(['region_index'], axis=1)
        cnv_regions = pd.concat([cnv_regions, regions], ignore_index=True)
    print("Join regions with copy-number variants of all chromosomes.")
    cnv_regions = cnv_regions[cnv_regions['cnv'] != 'haploid']
    cnv_regions = cnv_regions.round(2)
    cnv_regions['sample'] = sample_name

    print("Convert from 0-based to 1-based coordinates")    
    cnv_regions['start'] = cnv_regions['start'] + 1
    print("Save CNV regions to file.")
    output_path = Path(cnv_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cnv_regions.to_csv(output_path, sep='\t', index=False, header=True)
    
if __name__ == '__main__':
    cnv_calling()