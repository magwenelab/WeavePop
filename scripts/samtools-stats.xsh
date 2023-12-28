#/usr/bin/env xonsh -c

from pathlib import Path
import pandas as pd
import click

@click.command()
@click.argument("sample", type=str)
@click.argument("bamfile", type=str)
@click.argument("reference", type=str)
@click.argument("mapqfile", type=str)
@click.argument("covfile", type=str)

def stats(sample, bamfile, reference, mapqfile, covfile ): # Start definition of function with the sample as argument
    """This script runs 'samtools stats' on a sample for each chromosome, 
    extracts the MAPQ and COV sections for each chromosome and 
    combines the results of all chromosomes in 
    mapq.tsv and cov.tsv, respectively.
    
    SAMPLE is the sample name.
    BAMFILE is the path to the .bam
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
    out_cov = []
    for chromosome in chroms_list:
        mapq = $(samtools stats @(bamfile) @(chromosome) | grep ^MAPQ | cut -f 2-)
        mapq = pd.Series(list(mapq.split("\n")))
        mapq = chromosome + "\t" + mapq
        mapq = mapq.str.split("\t", expand = True)
        out_mapq.append(mapq)
        cov = $(samtools stats @(bamfile) @(chromosome) | grep ^COV | cut -f 2-)
        cov = pd.Series(list(cov.split("\n")))
        cov = chromosome + "\t" + cov
        cov= cov.str.split("\t", expand = True)
        out_cov.append(cov)
    quality = pd.concat(out_mapq)
    quality = quality.dropna()
    quality.columns = ["Chromosome", "MAPQ", "Count"]
    coverage = pd.concat(out_cov)
    coverage = coverage.dropna()
    coverage.columns = ["Chromosome", "Range", "Coverage", "Count"]
    quality.to_csv(mapqfile, index=False)
    coverage.to_csv(covfile, index=False)

if __name__ == "__main__":
    stats() # Runs the function stats