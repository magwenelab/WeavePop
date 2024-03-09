#/usr/bin/env xonsh -c

from pathlib import Path
import pandas as pd
import click
import os

@click.command()
@click.option("-s", "--sample", type=str, help="Sample name")
@click.option("-b", "--bamfile", type=click.Path(exists=True), help="Input BAM file")
@click.option("-g", "--bamgood", type=click.Path(exists=True), help="Good BAM file")
@click.option("-r", "--reference", type=click.Path(exists=True), help="Reference file")
@click.option("-cn", "--chrom_names", type=click.Path(), help="Chromosome names file")
@click.option("-m", "--mapq_out", type=click.Path(), help="Output MAPQ file")
@click.option("-c", "--cov_out", type=click.Path(), help="Output coverage file")
@click.option("-p", "--mapped_out", type=click.Path(), help="Output file with mapped reads metrics")

def stats(sample, bamfile, bamgood, reference, mapq_out, cov_out,mapped_out, chrom_names): # Start definition of function with the sample as argument
    chromosomes = pd.read_csv(chrom_names, names = ["Lineage", "Accession", "Chromosome"])
    lineage = os.path.splitext(Path(reference).name)[0]
    lin_chromosomes = chromosomes[chromosomes["Lineage"] == lineage]
    chroms_list = lin_chromosomes["Accession"].tolist()
    out_mapq = []
    out_cov_raw = []
    out_cov_good = []
    for chromosome in chroms_list:
        print("Analysing chromosome", chromosome)
        mapq = $(samtools stats @(bamfile) @(chromosome) | grep ^MAPQ | cut -f 2-)
        mapq = pd.Series(list(mapq.split("\n")))
        mapq = chromosome + "\t" + mapq
        mapq = mapq.str.split("\t", expand = True)
        out_mapq.append(mapq)
        cov_raw = $(samtools stats @(bamfile) @(chromosome) | grep ^COV | cut -f 2-)
        cov_raw = pd.Series(list(cov_raw.split("\n")))
        cov_raw = chromosome + "\t" + cov_raw
        cov_raw = cov_raw.str.split("\t", expand = True)
        out_cov_raw.append(cov_raw)
        cov_good = $(samtools stats @(bamgood) @(chromosome) | grep ^COV | cut -f 2-)
        cov_good = pd.Series(list(cov_good.split("\n")))
        cov_good = chromosome + "\t" + cov_good
        cov_good= cov_good.str.split("\t", expand = True)
        out_cov_good.append(cov_good)

    print("Saving MAPQ table")
    quality = pd.concat(out_mapq)
    quality = quality.dropna()
    quality.columns = ["Chromosome", "MAPQ", "Count"]
    quality.to_csv(mapq_out, index=False, sep = "\t")

    quality['MAPQ'] = quality['MAPQ'].astype(int)
    quality['Count'] = quality['Count'].astype(int)
    quality['MAPQ_20'] = pd.cut(quality['MAPQ'], bins=[-1, 19, 59, 100], labels=['MAPQ_20', 'MAPQ_20_59', 'MAPQ_60'])
    quality_sum = quality.groupby(['MAPQ_20']).agg({'Count': 'sum'}).reset_index()
    quality_sum.rename(columns={'Count': 'Count_bins'}, inplace=True)
    quality_sum['sample'] = sample
    quality_sum = quality_sum.round(2)
    quality_wider = quality_sum.pivot(index='sample', columns='MAPQ_20', values='Count_bins').reset_index()

    print("Saving coverage tables")
    coverage_raw = pd.concat(out_cov_raw)
    coverage_raw = coverage_raw.dropna()
    coverage_raw.columns = ["Chromosome", "Range", "Coverage", "Count_Raw"]
    coverage_raw.drop('Range', axis = 1, inplace = True)
    
    coverage_good = pd.concat(out_cov_good)
    coverage_good = coverage_good.dropna()
    coverage_good.columns = ["Chromosome", "Range", "Coverage", "Count_Good"]
    coverage_good.drop('Range', axis = 1, inplace = True)

    coverage = pd.merge(coverage_good, coverage_raw, on = ['Coverage', 'Chromosome'], how = 'outer')
    coverage.drop_duplicates( inplace=True)

    coverage.to_csv(cov_out, index=False, sep = "\t")

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
    stats_wider['percent_mapped'] = stats_wider['reads_mapped'] / stats_wider['raw_total_sequences'] * 100
    stats_wider['percent_paired'] = stats_wider['reads_properly_paired'] / stats_wider['raw_total_sequences'] * 100
    stats_wider = stats_wider.round(2)
    print("Joining mapped reads metrics with MAPQ metrics")
    stats_wider = pd.merge(stats_wider, quality_wider, on = 'sample', how = 'outer')
    print("Saving mapped reads metrics")
    stats_wider.to_csv(mapped_out, index=False, sep = "\t")
if __name__ == "__main__":
        stats() # Runs the function stats


