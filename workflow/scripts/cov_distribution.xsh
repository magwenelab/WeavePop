#/usr/bin/env xonsh -c

from pathlib import Path
import pandas as pd
import click
import os

@click.command()
@click.option("-s", "--sample", type=str, help="Sample name")
@click.option("-b", "--bamfile", type=click.Path(exists=True), help="Input BAM file")
@click.option("-g", "--bamgood", type=click.Path(exists=True), help="Good BAM file")
@click.option("-do", "--distribution_out", type=click.Path(), help="Output coverage file")
@click.option("-go", "--global_out", type=click.Path(), help="Output chromosomal coverage file")

def stats(sample, bamfile, bamgood, distribution_out, global_out): 
    chromosomes = $(samtools idxstats @(bamfile) | cut -f1 | grep -v "*")
    chromosomes = pd.Series(list(filter(None, chromosomes.split("\n"))))
    out_cov_raw = []
    out_cov_good = []
    print("Getting coverage distribution")
    for chromosome in chromosomes:
        print("Analysing chromosome", chromosome)
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

    print("Merging coverage distributions")
    coverage_raw = pd.concat(out_cov_raw)
    coverage_raw = coverage_raw.dropna()
    coverage_raw.columns = ["Accession", "Range", "Coverage", "Count_Raw"]
    coverage_raw.drop('Range', axis = 1, inplace = True)
    
    coverage_good = pd.concat(out_cov_good)
    coverage_good = coverage_good.dropna()
    coverage_good.columns = ["Accession", "Range", "Coverage", "Count_Good"]
    coverage_good.drop('Range', axis = 1, inplace = True)

    coverage = pd.merge(coverage_good, coverage_raw, on = ['Coverage', 'Accession'], how = 'outer')
    coverage.drop_duplicates( inplace=True)
    coverage.fillna(0, inplace=True)
    coverage.to_csv(distribution_out, index=False, sep = "\t")

    print("Get global coverage mode")
    coverage['Count_Good'] = coverage['Count_Good'].astype(int)
    global_cov = coverage.groupby('Coverage').agg({'Count_Good': 'sum'}).reset_index()
    print(global_cov)
    global_mode = global_cov.loc[global_cov['Count_Good'].idxmax()]['Coverage']
    print(global_mode)
    df = pd.DataFrame({'Sample': [sample], 'Global_Mode': [global_mode]})
    print(df)
    df.to_csv(global_out, index=False, sep = "\t")

if __name__ == "__main__":
        stats() 


