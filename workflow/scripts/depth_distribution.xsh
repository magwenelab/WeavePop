#/usr/bin/env xonsh -c

from pathlib import Path
import pandas as pd
import click
import os

@click.command()
@click.option("-s", "--sample", type=str, help="Sample name")
@click.option("-b", "--bamfile", type=click.Path(exists=True), help="Input BAM file")
@click.option("-g", "--bamgood", type=click.Path(exists=True), help="Imput good BAM file")
@click.option("-do", "--distribution_out", type=click.Path(), help="Distribution of depth output table")
@click.option("-go", "--chromosome_good_out", type=click.Path(), help="Global depth by chromosome of good alignments output table")
@click.option("-ro", "--chromosome_raw_out", type=click.Path(), help="Global depth by chromosome of raw alignments output table") 


def stats(sample, bamfile, bamgood, distribution_out, chromosome_good_out, chromosome_raw_out): 
    print("Getting chromosome names with samtools idxstats...")
    chromosomes = $(samtools idxstats @(bamfile) | cut -f1 | grep -v "*")
    chromosomes = pd.Series(list(filter(None, chromosomes.split("\n"))))
    out_depth_raw = []
    out_depth_good = []
    print("Getting depth distribution...")
    for chromosome in chromosomes:
        print("Analysing chromosome", chromosome, "...")
        depth_raw = $(samtools stats @(bamfile) @(chromosome) -c 1,1000000000,1| grep ^COV | cut -f 2-)
        depth_raw = pd.Series(list(depth_raw.split("\n")))
        depth_raw = chromosome + "\t" + depth_raw
        depth_raw = depth_raw.str.split("\t", expand = True)
        out_depth_raw.append(depth_raw)
        depth_good = $(samtools stats @(bamgood) @(chromosome) -c 1,1000000000,1| grep ^COV | cut -f 2-)
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

    print("Getting genome-wide depth...")
    distribution['count_raw'] = distribution['count_raw'].astype(int)
    distribution['count_good'] = distribution['count_good'].astype(int)
    distribution['depth'] = distribution['depth'].astype(int)
    global_cov = distribution.groupby('depth').agg({'count_good': 'sum'}).reset_index()
    global_mode = global_cov.loc[global_cov['count_good'].idxmax()]['depth']

    print("Getting mean and median depth by chromosome of good alignments...")
    def calculate_mean_median(group, count_quality):
        total_depth = (group['depth'] * group[count_quality]).sum()
        total_count = group[count_quality].sum()
        mean_depth = (total_depth / total_count).round(2)
        sorted_group = group.sort_values('depth')
        cutoff = sorted_group[count_quality].sum() / 2
        median_depth = sorted_group[sorted_group[count_quality].cumsum() >= cutoff].iloc[0]['depth']
        
        return pd.Series({'chrom_mean': mean_depth, 'chrom_median': median_depth})

    depth_by_chrom_good = distribution.groupby('accession').apply(calculate_mean_median, 'count_good', include_groups=False).reset_index()
    global_good = calculate_mean_median(distribution, 'count_good').reset_index()

    print("Gathering results of good alignments in dataframe ...")
    depth_by_chrom_good['global_mean'] = global_good[0][0]
    depth_by_chrom_good['global_median'] = global_good[0][1]
    depth_by_chrom_good['global_mode'] = global_mode

    print("Normalizing ...")
    genome_wide_depth = global_mode
    depth_by_chrom_good['norm_chrom_mean'] = depth_by_chrom_good['chrom_mean'] / genome_wide_depth
    depth_by_chrom_good['norm_chrom_median'] = depth_by_chrom_good['chrom_median'] / genome_wide_depth
    depth_by_chrom_good['norm_global_mean'] = depth_by_chrom_good['global_mean'] / genome_wide_depth
    depth_by_chrom_good['norm_global_median'] = depth_by_chrom_good['global_median'] / genome_wide_depth
    depth_by_chrom_good['sample'] = sample
    depth_by_chrom_good = depth_by_chrom_good.round(2)
    depth_by_chrom_good = depth_by_chrom_good[ ['sample'] + [ col for col in depth_by_chrom_good.columns if col != 'sample' ] ]

    print("Getting mean and median depth by chromosome for raw alignments...")
    depth_by_chrom_raw = distribution.groupby('accession').apply(calculate_mean_median, 'count_raw', include_groups=False).reset_index()
    global_raw = calculate_mean_median(distribution, 'count_raw').reset_index()

    print("Gathering results of raw alignments in dataframe...")
    depth_by_chrom_raw['global_mean'] = global_raw[0][0]
    depth_by_chrom_raw['global_median'] = global_raw[0][1]
    depth_by_chrom_raw['sample'] = sample
    depth_by_chrom_raw = depth_by_chrom_raw[ ['sample'] + [ col for col in depth_by_chrom_raw.columns if col != 'sample' ] ]

    print("Saving results...")
    depth_by_chrom_good.to_csv(chromosome_good_out, index=False, sep = "\t")
    depth_by_chrom_raw.to_csv(chromosome_raw_out, index=False, sep = "\t")

    print("Done!")
if __name__ == "__main__":
        stats() 




