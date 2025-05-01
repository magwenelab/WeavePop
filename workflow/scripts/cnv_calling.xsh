import pandas as pd
import io
import click
from scipy import ndimage
import numpy as np
from pathlib import Path
from itertools import product


@click.command()
@click.option('-di', '--depth_input', help='Path to BED file with depth of each window.', type=click.Path(exists=True))
@click.option('-ri', '--repeats_input', help='Path to BED file with coordinates of repetititve sequences of reference genome.', type=click.Path(exists=True))
@click.option('-ai', '--annotation_input', help='Path to TSV file with annotation of sample.', type=click.Path(exists=True))
@click.option('-ci', '--chromosome_input', help='Path to TSV file of chromosome lengths.', type=click.Path(exists=True))
@click.option('-sp', '--sample_name', help='Sample name as a string.', type=str)
@click.option('-wp', '--window_size', help='Size of windows in the depth BED file.', type=int)
@click.option('-dp', '--depth_threshold', help='Threshold to define copy number variation in smoothed normalzed depth.', type=click.types.FloatRange(min=0.0))
@click.option('-co', '--cnv_output', help='Path to output table of CNV calling.', type=click.Path())
@click.option('-mo', '--chromosome_metrics_output', help='Path to output table of CNV metrics per chromosome.', type=click.Path())
@click.option('-t', '--temp_dir', help='Path to temporary directory.', type=click.Path())
def cnv_calling(depth_input, repeats_input, annotation_input, chromosome_input, sample_name, window_size, depth_threshold, cnv_output, chromosome_metrics_output, temp_dir):

    print("Merging overlapping repeats and intersect with windows...")
    intersect = $(bedtools intersect -a @(depth_input) -b @(repeats_input) -wao)

    print("Reorganize intersection...")
    df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
    header = ['bed_accession', 'bed_start', 'bed_end','bed_depth', 'bed_norm_depth', 'bed_smooth_depth', 'r_accession', 'r_start', 'r_end', 'r_type', 'overlap_bp'] 
    df.columns = header

    print("Calculate overlap in base pairs...")
    df = df.drop(['r_accession', 'r_start', 'r_end'], axis=1)
    df['r_type_mix'] = df.groupby(['bed_accession', 'bed_start', 'bed_end', 'bed_depth', 'bed_norm_depth', 'bed_smooth_depth'])['r_type'].transform(lambda x: ','.join(x))
    df['overlap_bp_sum'] = df.groupby(['bed_accession', 'bed_start', 'bed_end', 'bed_depth', 'bed_norm_depth', 'bed_smooth_depth'])['overlap_bp'].transform('sum')
    df = df.drop(['r_type', 'overlap_bp'], axis=1)
    df = df.drop_duplicates()
    df.rename(columns={'r_type_mix': 'r_type', 'overlap_bp_sum': 'overlap_bp'}, inplace=True)    
    df = df.reset_index(drop=True)

    print("Calculating fraction of window with repetitive sequences...")
    repeats_fragments = df.copy()
    repeats_fragments['repeat_fraction'] = (repeats_fragments['overlap_bp'] / window_size).round(2)
    repeats_fragments.columns = repeats_fragments.columns.str.replace('bed_', '')
    repeats_fragments['sample'] = sample_name

    print("Defining copy-number of regions...")
    cnv_regions = pd.DataFrame()
    for accession in repeats_fragments['accession'].unique():
        regions_merged = repeats_fragments[repeats_fragments['accession'] == accession].copy()
        regions_merged.loc[:, 'cnv'] = 'single_copy'
        regions_merged = regions_merged.reset_index(drop=True)
        for i in range(len(regions_merged)):
            if regions_merged.loc[i, 'smooth_depth'] > 1 + depth_threshold:
                regions_merged.loc[i, 'cnv'] = "duplication"
            elif regions_merged.loc[i, 'smooth_depth'] < 1 - depth_threshold:
                regions_merged.loc[i, 'cnv'] = "deletion"
            else:
                regions_merged.loc[i, 'cnv'] = "single_copy"
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

    print("Rounding values and adding sample names...")
    cnv_regions = cnv_regions.round(2)
    cnv_regions['sample'] = sample_name

    print("Converting from 0-based to 1-based coordinates...")    
    cnv_regions['start'] = cnv_regions['start'] + 1

    print("Adding genetic features in CNV regions...")
    print("Filtering out single_copy regions from temporary dataframe...")
    temp_cnv_regions = cnv_regions[cnv_regions['cnv'] != 'single_copy']
    temp_cnv_file = Path(temp_dir) / f"{sample_name}_temp_cnv.bed"
    temp_cnv_regions.to_csv(temp_cnv_file, sep='\t', index=False, header=False)

    print("Intersecting with annotation...")
    annot_intersection = $(awk '$3 == "gene" {print $1,$4,$5,$10}' OFS='\t' @(annotation_input) | tr -d "'" | bedtools intersect -loj -a @(temp_cnv_file) -b stdin | bedtools merge -c 4,5,6,7,8,9,10,11,15 -o distinct)
    annot_intersection = pd.read_csv(io.StringIO(annot_intersection), sep='\t', header=None)
    temp_cnv_file.unlink()

    print("Naming and reording columns...")
    annot_header = ['accession', 'start', 'end', 'depth', 'norm_depth', 'smooth_depth', 'cnv', 'overlap_bp', 'region_size', 'repeat_fraction', 'sample', 'feature_id']
    annot_intersection.columns = annot_header
    annot_intersection['feature_id'] = annot_intersection['feature_id'].replace('.', np.nan)
    col_order = ['accession', 'start', 'end', 'cnv','region_size', 'depth', 'norm_depth', 'smooth_depth', 'repeat_fraction', 'overlap_bp', 'feature_id', 'sample']
    annot_intersection = annot_intersection[col_order]

    print("Saving CNV regions to file...")
    output_path = Path(cnv_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    annot_intersection.to_csv(output_path, sep='\t', index=False, header=True)


    print("Summarizing information for each chromosome and type of region...")
    regions = cnv_regions.groupby(['sample','accession', 'cnv']).agg({'norm_depth':['mean', 'median'],
                                                                        'smooth_depth':['mean', 'median'],   
                                                                        'region_size' :['sum', 'min', 'max', 'std'],     
                                                                        'start': ['size','min'],
                                                                        'end': 'max',
                                                                        },
                                                                        ).reset_index()

    regions.columns = ['sample','accession', 'cnv', 
                                        'norm_depth_mean',
                                        'norm_depth_median', 
                                        'smooth_depth_mean',
                                        'smooth_depth_median',
                                        'total_size_regions',
                                        'size_smallest_region',
                                        'size_largest_region',
                                        'std_regions_size',
                                        'n_regions',
                                        'first', 
                                        'last']

    print("Reading chromosome lengths...")
    chromosomes = pd.read_csv(chromosome_input, sep='\t', header=0)

    lineage = Path(repeats_input).parent.name

    chromosomes = chromosomes[chromosomes['lineage'] == lineage]


    print("Making sure that all types of CNV are present in all chromosomes...")
    all_combinations = pd.DataFrame(product([sample_name], 
                                            chromosomes['accession'].unique(), 
                                            ['duplication', 'deletion', 'single_copy']), 
                                    columns=['sample', 'accession', 'cnv'])

    all_combinations = pd.merge(all_combinations, 
                                    chromosomes, 
                                    how='left', 
                                    left_on='accession', 
                                    right_on='accession')


    print("Adding chromosome information to CNV regions...")
    regions_per_chromosome = pd.merge(all_combinations, regions, how='left', left_on=['sample','accession', 'cnv'], right_on=['sample','accession', 'cnv']) 
    regions_per_chromosome.loc[:, ['n_regions',  'total_size_regions']] = regions_per_chromosome.loc[:, ['n_regions', 'total_size_regions']].fillna(0).astype(int)

    print("Calculating coverage and span percentages...")
    regions_per_chromosome['coverage_percent'] = regions_per_chromosome['total_size_regions'] / regions_per_chromosome['length'] * 100
    regions_per_chromosome['span_percent'] = (regions_per_chromosome['last'] - regions_per_chromosome['first']) / regions_per_chromosome['length'] * 100
    regions_per_chromosome = regions_per_chromosome.drop(columns=['first', 'last'])
    regions_per_chromosome['span_percent'] = regions_per_chromosome['span_percent'].fillna(0)

    regions_per_chromosome = round(regions_per_chromosome, 2)

    print("Reorganizing columns...")
    regions_per_chromosome = regions_per_chromosome[['sample', 'lineage', 'accession', 'chromosome', 'length', 'cnv', 'n_regions', 'total_size_regions', 'coverage_percent', 'span_percent', 'size_smallest_region', 'size_largest_region', 'std_regions_size', 'norm_depth_mean', 'norm_depth_median', 'smooth_depth_mean', 'smooth_depth_median']]

    print("Saving summary of CNV regions to file...")
    regions_per_chromosome.to_csv(chromosome_metrics_output, sep='\t', index=False, header=True)

print("Done!")

if __name__ == '__main__':
    cnv_calling()





