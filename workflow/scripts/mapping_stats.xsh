#/usr/bin/env xonsh -c

from pathlib import Path
import pandas as pd
import click
import os

@click.command()
@click.option("-s", "--sample", type=str, help="Sample name")
@click.option("-b", "--bamfile", type=click.Path(exists=True), help="Input BAM file")
@click.option("-m", "--global_mode", type=click.Path(), help="Path to table with global mode of the sample")
@click.option("-l", "--low_mapq", type=int, help="Threshold of low MAPQ bin")
@click.option("-h", "--high_mapq", type=int, help="Threshold of high MAPQ bin")
@click.option("-pd", "--min_position_depth", type=int, help="Minimum depth of a position to be considered covered")
@click.option("-d", "--min_depth", type=int, help="Minimum percentage of genome-wide depth (Global mode)")
@click.option("-q", "--min_mapq", type=int, help="Minimum percentage of MAPQ high")
@click.option("-p", "--min_pp", type=int, help="Minimum percentage of properly paired reads")
@click.option("-c", "--min_coverage", type=int, help="Minimum percentage of genome covered")
@click.option("-o", "--output", type=click.Path(), help="Output file with mapped reads metrics")

def stats(sample, bamfile,  global_mode, low_mapq, high_mapq, min_position_depth, min_depth, min_mapq, min_pp, min_coverage, output):
    chromosomes = $(samtools idxstats @(bamfile) | cut -f1 | grep -v "*")
    chromosomes = pd.Series(list(filter(None, chromosomes.split("\n"))))

    print("Getting MAPQ distribution")
    out_mapq = []
    for chromosome in chromosomes:
        print("Analysing chromosome", chromosome)
        mapq = $(samtools stats @(bamfile) @(chromosome) | grep ^MAPQ | cut -f 2-)
        mapq = pd.Series(list(mapq.split("\n")))
        mapq = chromosome + "\t" + mapq
        mapq = mapq.str.split("\t", expand = True)
        out_mapq.append(mapq)

    print("Concatenating MAPQ distribution of all chromosomes")
    quality = pd.concat(out_mapq)
    quality = quality.dropna()
    quality.columns = ["Accession", "mapq", "count"]

    print("Binning MAPQ values")
    quality['mapq'] = quality['mapq'].astype(int)
    quality['count'] = quality['count'].astype(int)
    quality['low_mapq'] = pd.cut(quality['mapq'], bins=[-1, (low_mapq -1), (high_mapq -1), 100], labels=['low_mapq', 'intermediate_mapq', 'high_mapq'])
    quality_sum = quality.groupby(['low_mapq'], observed=False).agg({'count': 'sum'}).reset_index()
    quality_sum.rename(columns={'count': 'count_bins'}, inplace=True)
    quality_sum['sample'] = sample
    quality_sum = quality_sum.round(2)
    quality_wider = quality_sum.pivot(index='sample', columns='low_mapq', values='count_bins').reset_index()

    print("Getting mapped reads metrics")
    stats = $(samtools stats @(bamfile))
    sn_lines = [line.split('#')[0] for line in stats.split('\n') if line.startswith('SN')]
    sn_stats = '\n'.join(sn_lines)
    sn_stats = sn_stats.replace("SN", "")
    sn_stats = sn_stats.replace("\t", "")
    sn_stats = sn_stats.replace(":", "\t")
    sn_stats = sn_stats.split("\n")
    sn_stats = [x.split("\t") for x in sn_stats]
    sn_stats = pd.DataFrame(sn_stats)
    sn_stats.columns = ["Metric", "Value"]
    sn_stats['Value'] = sn_stats['Value'].astype(float)
    sn_stats['Metric'] = sn_stats['Metric'].str.replace(' ', '_')
    filtered_stats = sn_stats[sn_stats['Metric'].isin(['raw_total_sequences', 'reads_mapped', 'reads_properly_paired', 'average_quality'])].copy()
    filtered_stats['sample'] = sample
    stats_wider = filtered_stats.pivot(index='sample', columns='Metric', values='Value').reset_index()
    stats_wider['reads_unmapped'] = stats_wider['raw_total_sequences'] - stats_wider['reads_mapped']
    stats_wider['percent_unmapped'] = stats_wider['reads_unmapped'] / stats_wider['raw_total_sequences'] * 100
    stats_wider['percent_mapped'] = stats_wider['reads_mapped'] / stats_wider['raw_total_sequences'] * 100
    stats_wider['percent_paired'] = stats_wider['reads_properly_paired'] / stats_wider['raw_total_sequences'] * 100
    stats_wider['reads_only_mapped'] = stats_wider['reads_mapped'] - stats_wider['reads_properly_paired']
    stats_wider['percent_only_mapped'] = stats_wider['reads_only_mapped'] / stats_wider['raw_total_sequences'] * 100
    stats_wider['percent_properly_paired'] = stats_wider['reads_properly_paired'] / stats_wider['raw_total_sequences'] * 100

    print("Joining mapped reads metrics with MAPQ metrics")
    stats_wider = pd.merge(stats_wider, quality_wider, on = 'sample', how = 'outer')
    stats_wider['percent_low_mapq'] = stats_wider['low_mapq'] / stats_wider['reads_mapped'] * 100
    stats_wider['percent_inter_mapq'] = stats_wider['intermediate_mapq'] / stats_wider['reads_mapped'] * 100
    stats_wider['percent_high_mapq'] = stats_wider['high_mapq'] / stats_wider['reads_mapped'] * 100
    stats_wider = stats_wider.round(2)

    print("Joining mapped reads metrics with global mode")
    global_mode = pd.read_csv(global_mode, sep = "\t", header = 0)
    global_mode = global_mode['Global_Mode'][0]
    stats_wider['genome-wide_depth'] = global_mode
    
    print("Calculating coverage")
    depth = $(samtools depth -a @(bamfile))
    depth_string_list = depth.split('\n')
    depth_df = pd.DataFrame([x.split('\t') for x in depth_string_list], columns=['chrom', 'pos', 'depth'])
    depth_df = depth_df.iloc[:-1]
    depth_df['depth'] = depth_df['depth'].astype(int)
    coverage = (depth_df[depth_df['depth'] > min_position_depth].shape[0]/ depth_df.shape[0]) * 100
    stats_wider['percent_covered'] = round(coverage, 2)

    print("Adding quality warning flag")
    stats_wider['mapq_warning'] = stats_wider.apply(lambda row: "MAPQ-Low" if row['percent_high_mapq'] < min_mapq else None, axis=1)
    stats_wider['pp_warning'] = stats_wider.apply(lambda row: "Properly-paired-Low" if row['percent_properly_paired'] < min_pp else None, axis=1)
    stats_wider['depth_warning'] = stats_wider.apply(lambda row: "Depth-Low" if row['genome-wide_depth'] <= min_depth else None, axis=1)
    stats_wider['coverage_warning'] = stats_wider.apply(lambda row: "Coverage-Low" if row['percent_covered'] < min_coverage else None, axis=1)

    print("Joining warnings")
    stats_wider['quality_warning'] = stats_wider[['mapq_warning', 'pp_warning', 'depth_warning','coverage_warning']].apply(lambda x: '_'.join(x.dropna()), axis=1)
    stats_wider = stats_wider.drop(columns = ['mapq_warning', 'pp_warning', 'depth_warning','coverage_warning'])
    
    print("Saving mapped reads metrics")
    stats_wider.to_csv(output, index=False, sep = "\t")
if __name__ == "__main__":
        stats() 


