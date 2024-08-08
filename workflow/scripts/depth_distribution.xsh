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
    chromosomes = $(samtools idxstats @(bamfile) | cut -f1 | grep -v "*")
    chromosomes = pd.Series(list(filter(None, chromosomes.split("\n"))))
    out_depth_raw = []
    out_depth_good = []
    print("Getting depth distribution")
    for chromosome in chromosomes:
        print("Analysing chromosome", chromosome)
        depth_raw = $(samtools stats @(bamfile) @(chromosome) | grep ^COV | cut -f 2-)
        depth_raw = pd.Series(list(depth_raw.split("\n")))
        depth_raw = chromosome + "\t" + depth_raw
        depth_raw = depth_raw.str.split("\t", expand = True)
        out_depth_raw.append(depth_raw)
        depth_good = $(samtools stats @(bamgood) @(chromosome) | grep ^COV | cut -f 2-)
        depth_good = pd.Series(list(depth_good.split("\n")))
        depth_good = chromosome + "\t" + depth_good
        depth_good= depth_good.str.split("\t", expand = True)
        out_depth_good.append(depth_good)

    print("Merging depth distributions")
    distribution_raw = pd.concat(out_depth_raw)
    distribution_raw = distribution_raw.dropna()
    distribution_raw.columns = ["Accession", "Range", "Depth", "Count_Raw"]
    distribution_raw.drop('Range', axis = 1, inplace = True)
    
    distribution_good = pd.concat(out_depth_good)
    distribution_good = distribution_good.dropna()
    distribution_good.columns = ["Accession", "Range", "Depth", "Count_Good"]
    distribution_good.drop('Range', axis = 1, inplace = True)

    distribution = pd.merge(distribution_good, distribution_raw, on = ['Depth', 'Accession'], how = 'outer')
    distribution.drop_duplicates( inplace=True)
    distribution.fillna(0, inplace=True)
    distribution.to_csv(distribution_out, index=False, sep = "\t")

    print("Get global depth mode")
    distribution['Count_Raw'] = distribution['Count_Raw'].astype(int)
    distribution['Count_Good'] = distribution['Count_Good'].astype(int)
    distribution['Depth'] = distribution['Depth'].astype(int)
    global_cov = distribution.groupby('Depth').agg({'Count_Good': 'sum'}).reset_index()
    global_mode = global_cov.loc[global_cov['Count_Good'].idxmax()]['Depth']

    print("Get mean and median depth by chromosome")
    def calculate_mean_median(group, quality):
        total_depth = (group['Depth'] * group[quality]).sum()
        total_count = group[quality].sum()
        mean_depth = (total_depth / total_count).round(2)
        sorted_group = group.sort_values('Depth')
        cutoff = sorted_group[quality].sum() / 2
        median_depth = sorted_group[sorted_group[quality].cumsum() >= cutoff].iloc[0]['Depth']
        
        return pd.Series({'Chrom_Mean': mean_depth, 'Chrom_Median': median_depth})

    depth_by_chrom_good = distribution.groupby('Accession').apply(calculate_mean_median, 'Count_Good', include_groups=False).reset_index()
    global_good = calculate_mean_median(distribution, 'Count_Good').reset_index()

    depth_by_chrom_good['Global_Mean'] = global_good[0][0]
    depth_by_chrom_good['Global_Median'] = global_good[0][1]
    depth_by_chrom_good['Global_Mode'] = global_mode
    depth_by_chrom_good['Norm_Chrom_Mean'] = depth_by_chrom_good['Chrom_Mean'] / global_mode
    depth_by_chrom_good['Norm_Chrom_Median'] = depth_by_chrom_good['Chrom_Median'] / global_mode
    depth_by_chrom_good['Norm_Global_Mean'] = depth_by_chrom_good['Global_Mean'] / global_mode
    depth_by_chrom_good['Norm_Global_Median'] = depth_by_chrom_good['Global_Median'] / global_mode
    depth_by_chrom_good['Sample'] = sample
    depth_by_chrom_good = depth_by_chrom_good.round(2)
    depth_by_chrom_good = depth_by_chrom_good[ ['Sample'] + [ col for col in depth_by_chrom_good.columns if col != 'Sample' ] ]


    depth_by_chrom_raw = distribution.groupby('Accession').apply(calculate_mean_median, 'Count_Raw', include_groups=False).reset_index()
    global_raw = calculate_mean_median(distribution, 'Count_Raw').reset_index()

    depth_by_chrom_raw['Global_Mean'] = global_raw[0][0]
    depth_by_chrom_raw['Global_Median'] = global_raw[0][1]
    depth_by_chrom_raw['Sample'] = sample
    depth_by_chrom_raw = depth_by_chrom_raw[ ['Sample'] + [ col for col in depth_by_chrom_raw.columns if col != 'Sample' ] ]


    depth_by_chrom_good.to_csv(chromosome_good_out, index=False, sep = "\t")
    depth_by_chrom_raw.to_csv(chromosome_raw_out, index=False, sep = "\t")

if __name__ == "__main__":
        stats() 




