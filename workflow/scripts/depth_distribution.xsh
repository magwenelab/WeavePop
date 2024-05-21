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
@click.option("-go", "--global_out", type=click.Path(), help="Global depth mode output table")

def stats(sample, bamfile, bamgood, distribution_out, global_out): 
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
    distribution['Count_Good'] = distribution['Count_Good'].astype(int)
    global_cov = distribution.groupby('Depth').agg({'Count_Good': 'sum'}).reset_index()
    print(global_cov)
    global_mode = global_cov.loc[global_cov['Count_Good'].idxmax()]['Depth']
    print(global_mode)
    df = pd.DataFrame({'Sample': [sample], 'Global_Mode': [global_mode]})
    print(df)
    df.to_csv(global_out, index=False, sep = "\t")

if __name__ == "__main__":
        stats() 


