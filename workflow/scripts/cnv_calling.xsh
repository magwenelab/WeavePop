import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import io
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import product
from scipy import ndimage

depth_input=snakemake.input.depth
repeats_input=snakemake.input.repeats
annotation_input=snakemake.input.annotation
chromosome_input=snakemake.input.chrom_length

depth_output=snakemake.output.windows
cnv_output=snakemake.output.cnvs
chromosome_metrics_output=snakemake.output.chrom_cnvs

sample_name=snakemake.wildcards.sample
smoothing_size=snakemake.params.smoothing_size
window_size=snakemake.params.window_size
depth_threshold=snakemake.params.depth_threshold

temp_dir=snakemake.resources.tmpdir


print("Merging overlapping repeats and intersect with windows...")
intersect = $(bedtools intersect -a @(depth_input) -b @(repeats_input) -wao 2>> @(log_file))

print("Reorganize intersection...")
df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)
header = ['bed_accession', 'bed_start', 'bed_end','bed_depth', 'r_accession', 'r_start', 'r_end', 'repeat_types', 'repeat_overlap_bp'] 
df.columns = header

print("Calculate overlap in base pairs...")
df = df.drop(['r_accession', 'r_start', 'r_end'], axis=1)
df['repeat_types_mix'] = df.groupby(['bed_accession', 'bed_start', 'bed_end', 'bed_depth'])['repeat_types'].transform(lambda x: ','.join(x.str.split(',').explode().unique()))
df['repeat_overlap_bp_sum'] = df.groupby(['bed_accession', 'bed_start', 'bed_end', 'bed_depth'])['repeat_overlap_bp'].transform('sum')
df = df.drop(['repeat_types', 'repeat_overlap_bp'], axis=1)
df = df.drop_duplicates()
df.rename(columns={'repeat_types_mix': 'repeat_types', 'repeat_overlap_bp_sum': 'repeat_overlap_bp'}, inplace=True)    
df = df.reset_index(drop=True)
df.columns = df.columns.str.replace('bed_', '')

print("Normalizing depth by median depth...")
genome_median_depth = df['depth'].median()
df.loc[:,'norm_depth'] = df['depth'] / genome_median_depth
df.loc[:,'norm_depth'] = df.loc[:,'norm_depth']

print("Smoothing depth values...")
cov_array = np.array(df["norm_depth"])
smoothed_array = ndimage.median_filter(cov_array, size=smoothing_size)
df.loc[:,'smooth_depth']=pd.Series(smoothed_array)

print("Reorganizing columns...")
df = df[['accession', 'start', 'end', 'depth', 'norm_depth', 'smooth_depth', 'repeat_types', 'repeat_overlap_bp']]
df = df.round(2)
df.to_csv(depth_output, sep='\t', index=False, header=True)

print("Defining copy-number of windows...")
cnv_windows = pd.DataFrame()

for accession in df['accession'].unique():
    windows = df[df['accession'] == accession].copy()
    windows.loc[:, 'cnv'] = 'single_copy'
    windows = windows.reset_index(drop=True)

    for i in range(len(windows)):
        if windows.loc[i, 'smooth_depth'] > 1 + depth_threshold:
            windows.loc[i, 'cnv'] = "duplication"
        elif windows.loc[i, 'smooth_depth'] < 1 - depth_threshold:
            windows.loc[i, 'cnv'] = "deletion"
        else:
            windows.loc[i, 'cnv'] = "single_copy"

    windows.loc[:,'region_index'] = 1

    for i in range(1, len(windows)):
        if windows.loc[i, 'cnv'] == windows.loc[i - 1, 'cnv']:
            windows.loc[i, 'region_index'] = windows.loc[i - 1, 'region_index']
        else:
            windows.loc[i, 'region_index'] = windows.loc[i - 1, 'region_index'] + 1

    cnv_windows = pd.concat([cnv_windows, windows], ignore_index=True)

print("Merging windows with same CNV into regions...")
cnv_regions = cnv_windows.groupby(['accession','region_index']).agg(start=('start', 'first'), 
                                                                    end=('end', 'last'),
                                                                    depth = ('depth', 'median'), 
                                                                    norm_depth=('norm_depth', 'median'), 
                                                                    smooth_depth=('smooth_depth', 'median'),
                                                                    cnv=('cnv', 'first'), 
                                                                    repeat_overlap_bp=('repeat_overlap_bp', 'sum')).reset_index()

cnv_regions['region_size'] = cnv_regions['end'] - cnv_regions['start']
cnv_regions['repeat_fraction'] = (cnv_regions['repeat_overlap_bp'] / cnv_regions['region_size']).round(2)
cnv_regions = cnv_regions.drop(['region_index'], axis=1)

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
annot_header = ['accession', 'start', 'end', 'depth', 'norm_depth', 'smooth_depth', 'cnv', 'repeat_overlap_bp', 'region_size', 'repeat_fraction', 'sample', 'feature_id']
annot_intersection.columns = annot_header
annot_intersection['feature_id'] = annot_intersection['feature_id'].replace('.', np.nan)
col_order = ['sample', 'accession', 'start', 'end', 'cnv','region_size', 'depth', 'norm_depth', 'smooth_depth', 'repeat_fraction', 'repeat_overlap_bp', 'feature_id']
annot_intersection = annot_intersection[col_order]

print("Saving CNV regions to file...")
output_path = Path(cnv_output)
output_path.parent.mkdir(parents=True, exist_ok=True)
annot_intersection.to_csv(output_path, sep='\t', index=False, header=True)

print("Calculating depth of chromosomes...")
chromosome_depth = cnv_windows.groupby('accession')['depth'].median().round(2).reset_index()
chromosome_depth.columns = ['accession', 'chrom_median']

print("Normalizing chromosome depth...")
genome_median_depth = cnv_windows['depth'].median().round(4)

chromosome_depth['genome_median_depth'] = genome_median_depth
chromosome_depth['norm_chrom_median'] = chromosome_depth['chrom_median'] / genome_median_depth
chromosome_depth['norm_chrom_median'] = chromosome_depth['norm_chrom_median'].round(2)

print("Summarizing information for each chromosome and type of region...")
summary_windows = cnv_windows.groupby(['accession', 'cnv']).agg({'norm_depth':['mean', 'median'],
                                                        'smooth_depth':['mean', 'median'],   
                                                        }).reset_index()
summary_windows.columns = ['accession', 'cnv', 
                                    'norm_depth_mean',
                                    'norm_depth_median', 
                                    'smooth_depth_mean',
                                    'smooth_depth_median']

summary_regions = cnv_regions.groupby(['accession', 'cnv']).agg({'region_size' :['sum', 'min', 'max', 'std'],     
                                                        'start': ['size','min'],
                                                        'end': 'max',
                                                        }).reset_index()
summary_regions.columns = ['accession', 'cnv', 
                                    'total_size_regions',
                                    'size_smallest_region',
                                    'size_largest_region',
                                    'std_regions_size',
                                    'n_regions',
                                    'first', 
                                    'last']

summary = pd.merge(summary_windows, summary_regions, how='left', left_on=['accession', 'cnv'], right_on=['accession', 'cnv'])
summary['sample'] = sample_name

print("Reading chromosome lengths...")
chromosomes = pd.read_csv(chromosome_input, sep=',', header=0)

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
regions_per_chromosome = pd.merge(all_combinations, summary, 
                                how='left', 
                                left_on=['sample','accession', 'cnv'], 
                                right_on=['sample','accession', 'cnv']) 
regions_per_chromosome.loc[:, ['n_regions',  'total_size_regions']] = regions_per_chromosome.loc[:, ['n_regions', 'total_size_regions']].fillna(0).astype(int)

print("Calculating coverage and span percentages...")
regions_per_chromosome['coverage_percent'] = regions_per_chromosome['total_size_regions'] / regions_per_chromosome['length'] * 100
regions_per_chromosome['span_percent'] = (regions_per_chromosome['last'] - regions_per_chromosome['first']) / regions_per_chromosome['length'] * 100
regions_per_chromosome = regions_per_chromosome.drop(columns=['first', 'last'])
regions_per_chromosome['span_percent'] = regions_per_chromosome['span_percent'].fillna(0)

regions_per_chromosome = round(regions_per_chromosome, 2)

print("Reorganizing columns...")
regions_per_chromosome = regions_per_chromosome[['sample', 'cnv', 'accession', 'n_regions', 'total_size_regions', 'coverage_percent', 'span_percent', 'size_smallest_region', 'size_largest_region', 'std_regions_size', 'norm_depth_mean', 'norm_depth_median', 'smooth_depth_mean', 'smooth_depth_median']]

print("Adding chromosome depth to summary...")
regions_per_chromosome = pd.merge(regions_per_chromosome, chromosome_depth, how='left', left_on='accession', right_on='accession')

print("Saving summary of CNV regions to file...")
regions_per_chromosome.to_csv(chromosome_metrics_output, sep='\t', index=False, header=True)

print("Done!")





