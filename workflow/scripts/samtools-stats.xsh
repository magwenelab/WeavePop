#/usr/bin/env xonsh -c

from pathlib import Path
import pandas as pd
import click

@click.command()
@click.argument("sample", type=str)
@click.argument("bamfile", type=str)
@click.argument("bamgood", type=str)
@click.argument("reference", type=str)
@click.argument("mapqfile", type=str)
@click.argument("covfile", type=str)

def stats(sample, bamfile, bamgood, reference, mapqfile, covfile ): # Start definition of function with the sample as argument
    """This script runs 'samtools stats' on a sample for each chromosome, 
    extracts the MAPQ and COV sections for each chromosome and 
    combines the results of all chromosomes in 
    mapq.tsv and cov.tsv, respectively.
    
    SAMPLE is the sample name.
    BAMFILE is the path to the .bam
    BAMGOOD is the path to the bam file of only good quality alignments
    REFERENCE is the path to the ref.fa
    MAPQFILE is the path to the output table with the MAPQ results
    COVFILE is the path to the output table with the MAPQ results
    """
    chroms = $(grep chromosome @(reference))
    chroms_list = chroms.split("\n")
    chroms_list = list(filter(None, chroms_list))
    chroms_list = [i.split(' ',1)[0] for i in chroms_list]
    chroms_list = [i.replace('>', '') for i in chroms_list]
    out_mapq = []
    out_cov_raw = []
    out_cov_good = []
    for chromosome in chroms_list:
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

    quality = pd.concat(out_mapq)
    quality = quality.dropna()
    quality.columns = ["Chromosome", "MAPQ", "Count"]
    quality.to_csv(mapqfile, index=False)

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

    coverage.to_csv(covfile, index=False)

if __name__ == "__main__":
    stats() # Runs the function stats

