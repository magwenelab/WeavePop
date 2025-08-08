import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import os
from pathlib import Path
import pandas as pd

sample=snakemake.wildcards.unf_sample

bamfile=snakemake.input.bam
bamgood=snakemake.input.bam_good

distribution_out=snakemake.output.distrib
depth_summary_out=snakemake.output.summary

print("Getting chromosome names with samtools idxstats...")
chromosomes = $(samtools idxstats @(bamfile) | cut -f1 | grep -v "*" 2>> @(log_file))
chromosomes = pd.Series(list(filter(None, chromosomes.split("\n"))))
out_depth_raw = []
out_depth_good = []
print("Getting depth distribution...")
for chromosome in chromosomes:
    print("Analysing chromosome", chromosome, "...")
    depth_raw = $(samtools stats @(bamfile) @(chromosome) -c 1,1000000000,1| grep ^COV | cut -f 2- 2>> @(log_file))
    depth_raw = pd.Series(list(depth_raw.split("\n")))
    depth_raw = chromosome + "\t" + depth_raw
    depth_raw = depth_raw.str.split("\t", expand = True)
    out_depth_raw.append(depth_raw)
    depth_good = $(samtools stats @(bamgood) @(chromosome) -c 1,1000000000,1| grep ^COV | cut -f 2- 2>> @(log_file))
    depth_good = pd.Series(list(depth_good.split("\n")))
    depth_good = chromosome + "\t" + depth_good
    depth_good= depth_good.str.split("\t", expand = True)
    out_depth_good.append(depth_good)

print("Merging depth distributions...")
distribution_raw = pd.concat(out_depth_raw)
distribution_raw = distribution_raw.dropna()
distribution_raw.columns = ["accession", "range", "depth", "count_raw"]
distribution_raw.drop('range', axis = 1, inplace = True)

distribution_good = pd.concat(out_depth_good)
distribution_good = distribution_good.dropna()
distribution_good.columns = ["accession", "range", "depth", "count_good"]
distribution_good.drop('range', axis = 1, inplace = True)

distribution = pd.merge(distribution_good, distribution_raw, on = ['depth', 'accession'], how = 'outer')
distribution.drop_duplicates( inplace=True)
distribution.fillna(0, inplace=True)
distribution.to_csv(distribution_out, index=False, sep = "\t")

distribution['count_raw'] = distribution['count_raw'].astype(int)
distribution['count_good'] = distribution['count_good'].astype(int)
distribution['depth'] = distribution['depth'].astype(int)
distribution['sample'] = sample

print("Getting mean and median depth by chromosome and genome-wide...")
def calculate_mean_median(group, count_quality):
    total_depth = (group['depth'] * group[count_quality]).sum()
    total_count = group[count_quality].sum()
    mean_depth = (total_depth / total_count).round(2)
    sorted_group = group.sort_values('depth')
    cutoff = sorted_group[count_quality].sum() / 2
    median_depth = sorted_group[sorted_group[count_quality].cumsum() >= cutoff].iloc[0]['depth']
    
    return pd.Series({'mean': mean_depth, 'median': median_depth})

depth_by_chrom_good = distribution.groupby('accession').apply(calculate_mean_median, 'count_good', include_groups=False).reset_index()
depth_by_chrom_raw = distribution.groupby('accession').apply(calculate_mean_median, 'count_raw', include_groups=False).reset_index()
depth_genome_good = distribution.groupby('sample').apply(calculate_mean_median, 'count_good', include_groups=False).reset_index()
depth_genome_raw = distribution.groupby('sample').apply(calculate_mean_median, 'count_raw', include_groups=False).reset_index()

print("Gathering results of good alignments in dataframe ...")

summary_depth = pd.merge(depth_by_chrom_good, depth_by_chrom_raw, on = ['accession'], how = 'outer')
summary_depth.columns = ['accession', 'chrom_mean_good', 'chrom_median_good', 'chrom_mean_raw', 'chrom_median_raw']
summary_depth['genome_mean_good'] = depth_genome_good['mean'][0]
summary_depth['genome_median_good'] = depth_genome_good['median'][0]
summary_depth['genome_mean_raw'] = depth_genome_raw['mean'][0]
summary_depth['genome_median_raw'] = depth_genome_raw['median'][0]
summary_depth['sample'] = sample
summary_depth = summary_depth[ ['sample'] + [ col for col in summary_depth.columns if col != 'sample' ] ]

print("Saving summary table...")
summary_depth.to_csv(depth_summary_out, index=False, sep = "\t")

print("Done!")




