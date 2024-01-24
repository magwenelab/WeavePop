#!/usr/bin/env python

from collections import defaultdict
from pathlib import Path
import sys
import pandas as pd
import click


# Create dataframe from SAMPLESTABLE and REFGENOMETABLE with the sample, read filenames, lineage and reference assembly file information
@click.command()
@click.option('-s', '--samplestable', 'samplestable', 
              required=True, type=str, 
              help = "Sample metadata table with columns sample and group")
@click.option('-o', '--output', 'output', 
              required=True, type=str, 
              help = "File name of output table with the samples and the corresponding reference genome files.")
@click.option('-f1', '--suffix1', 'suffix1', 
              required=True, type=str, 
              help = "Suffix of the forward fastq filenames.")
@click.option('-f2', '--suffix2', 'suffix2', 
              required=True, type=str, 
              help = "Suffix of the referse fastq filenames.")
def getreference(samplestable,output, suffix1, suffix2):
    sample_data = pd.read_csv(samplestable)
    sample_data = sample_data[["sample","lineage"]]

    d={'lineage': sample_data["lineage"],
       'sample': sample_data["sample"], 
       'fastq1': sample_data["sample"]+ suffix1, 
       'fastq2': sample_data["sample"]+ suffix2, 
       'refgenome': sample_data["lineage"]+ ".fasta"}
    
    df = pd.DataFrame(data=d)

    filepath = Path(output)  
    df.to_csv(filepath, index=False)

if __name__ == "__main__":
    getreference() 