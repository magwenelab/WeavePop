import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import os
from pathlib import Path
import pandas as pd

sample=snakemake.wildcards.unf_sample
bamfile=snakemake.input.bam
depth=snakemake.input.depth
raw=snakemake.input.raw
filt=snakemake.input.filt
low_mapq=snakemake.params.low_mapq
high_mapq=snakemake.params.high_mapq
min_mapq=snakemake.params.min_mapq

output=snakemake.output[0]

print("Input parameters:")
print(f"Sample: {sample}")
print(f"BAM file: {bamfile}")
print(f"Depth file: {depth}")
print(f"Raw VCF file: {raw}")
print(f"Filtered VCF file: {filt}")
print(f"Low MAPQ threshold: {low_mapq}")
print(f"High MAPQ threshold: {high_mapq}")
print(f"Minimum MAPQ: {min_mapq}")
print(f"Output file: {output}")

print("Getting chromosome names...")
chromosomes = $(samtools idxstats @(bamfile) | cut -f1 | grep -v "*" 2>> @(log_file))
chromosomes = pd.Series(list(filter(None, chromosomes.split("\n"))))

print("Getting MAPQ distribution...")
quality = pd.DataFrame(columns=["accession", "mapq", "count"])
for chromosome in chromosomes:
    print("Analysing chromosome", chromosome, "...")
    mapq = $(samtools stats @(bamfile) @(chromosome) | grep ^MAPQ | cut -f 2- 2>> @(log_file))
    mapq = mapq.split("\n")
    mapq = [x.split("\t") for x in mapq if x]
    mapq_df = pd.DataFrame(mapq, columns=["mapq", "count"])
    mapq_df["accession"] = chromosome
    quality = pd.concat([quality, mapq_df], ignore_index=True)
    del mapq, mapq_df

print("Binning MAPQ values...")
quality['mapq'] = quality['mapq'].astype(int)
quality['count'] = quality['count'].astype(int)
quality['low_mapq'] = pd.cut(quality['mapq'], bins=[-1, (low_mapq -1), (high_mapq -1), 100], labels=['low_mapq', 'intermediate_mapq', 'high_mapq'])
quality_sum = quality.groupby(['low_mapq'], observed=False).agg({'count': 'sum'}).reset_index()
quality_sum.rename(columns={'count': 'count_bins'}, inplace=True)
quality_sum['sample'] = sample
quality_sum = quality_sum.round(2)
quality_sum = quality_sum.pivot(index='sample', columns='low_mapq', values='count_bins').reset_index()

print("Getting mapped reads metrics...")
stats = $(samtools stats @(bamfile) -c 1,1000000000,1 2>> @(log_file))
stats = [line.split('#')[0] for line in stats.split('\n') if line.startswith('SN')]
stats = '\n'.join(stats)
stats = stats.replace("SN", "")
stats = stats.replace("\t", "")
stats = stats.replace(":", "\t")
stats = stats.split("\n")
stats = [x.split("\t") for x in stats]
stats = pd.DataFrame(stats)
stats.columns = ["metric", "value"]
stats['value'] = stats['value'].astype(float)
stats['metric'] = stats['metric'].str.replace(' ', '_')
stats = stats[stats['metric'].isin(['raw_total_sequences', 'reads_mapped', 'reads_properly_paired', 'average_quality'])].copy()
stats['sample'] = sample
stats = stats.pivot(index='sample', columns='metric', values='value').reset_index()
stats['reads_unmapped'] = stats['raw_total_sequences'] - stats['reads_mapped']
stats['percent_unmapped'] = stats['reads_unmapped'] / stats['raw_total_sequences'] * 100
stats['percent_mapped'] = stats['reads_mapped'] / stats['raw_total_sequences'] * 100
stats['percent_paired'] = stats['reads_properly_paired'] / stats['raw_total_sequences'] * 100
stats['reads_only_mapped'] = stats['reads_mapped'] - stats['reads_properly_paired']
stats['percent_only_mapped'] = stats['reads_only_mapped'] / stats['raw_total_sequences'] * 100
stats['percent_properly_paired'] = stats['reads_properly_paired'] / stats['raw_total_sequences'] * 100

print("Joining mapped reads metrics with MAPQ metrics...")
stats = pd.merge(stats, quality_sum, on = 'sample', how = 'outer')
stats['percent_low_mapq'] = stats['low_mapq'] / stats['reads_mapped'] * 100
stats['percent_inter_mapq'] = stats['intermediate_mapq'] / stats['reads_mapped'] * 100
stats['percent_high_mapq'] = stats['high_mapq'] / stats['reads_mapped'] * 100
stats = stats.round(2)

print("Joining mapped reads metrics with genome-wide depth...")
depth_stats = pd.read_csv(depth, sep = "\t", header = 0)
depth_stats = depth_stats [['sample','genome_mean_good',  'genome_median_good', 'genome_mean_raw', 'genome_median_raw']]
depth_stats.columns = ['sample', 'genome_mean_depth_good', 'genome_median_depth_good', 'genome_mean_depth_raw', 'genome_median_depth_raw']
depth_stats.drop_duplicates(inplace=True)

stats = pd.merge(stats, depth_stats , on = 'sample', how = 'outer')

print("Calculating coverage...")
coverage_good = $(samtools coverage @(bamfile) --min-MQ @(min_mapq) 2>> @(log_file))
coverage_good = coverage_good.split('\n')
coverage_good = [x.split('\t') for x in coverage_good if x]
coverage_good = pd.DataFrame(coverage_good[1:], columns=coverage_good[0])
coverage_good = (coverage_good['covbases'].astype(int).sum() / coverage_good['endpos'].astype(int).sum()) * 100
stats['coverage_good'] = round(coverage_good, 2)

coverage_raw = $(samtools coverage @(bamfile) 2>> @(log_file))
coverage_raw = coverage_raw.split('\n')
coverage_raw = [x.split('\t') for x in coverage_raw if x]
coverage_raw = pd.DataFrame(coverage_raw[1:], columns=coverage_raw[0])
coverage_raw = (coverage_raw['covbases'].astype(int).sum() / coverage_raw['endpos'].astype(int).sum()) * 100
stats['coverage_raw'] = round(coverage_raw, 2)

print("Obtaining number of raw and filtered variants...")

n_raw = int($(grep -v "#" @(raw) | wc -l))
n_filtered = int($(grep -v "#" @(filt) | wc -l))
n_removed = n_raw - n_filtered
percent_filtered = n_filtered / n_raw * 100
percent_removed = n_removed / n_raw * 100

stats['n_raw_vars'] = n_raw
stats['n_filtered_vars'] = n_filtered
stats['n_removed_vars'] = n_removed
stats['percent_filtered_vars'] = round(percent_filtered,2)
stats['percent_removed_vars'] = round(percent_removed,2)

print("Saving mapping stats table...")
stats.to_csv(output, index=False, sep = "\t")

print("Done!")

