from collections import defaultdict
from pathlib import Path
import sys
import pandas as pd
import click


# Create dataframe from SAMPLESTABLE and REFGENOMETABLE with the sample, read filenames, lineage and reference assembly file information
@click.command()
@click.option('-s', '--samplestable', 'samplestable', required=True, type=str, help = "Sample metadata table with columns sample and group")
@click.option('-o', '--output', 'output', required=True, type=str, help = "File name of output table with the samples and the corresponding reference genome files.")
@click.option('-f1', '--suffix1', 'suffix1', required=True, type=str, help = "Suffix of the forward fastq filenames.")
@click.option('-f2', '--suffix2', 'suffix2', required=True, type=str, help = "Suffix of the referse fastq filenames.")

def getreference(samplestable,output, suffix1, suffix2):
    lineage = pd.read_csv(samplestable)
    lineage = lineage[["sample","group"]]

    d={'group': lineage["group"],'sample': lineage["sample"], 'file1': lineage["sample"]+ suffix1, 'file2': lineage["sample"]+ suffix2, 'refgenome' : lineage["group"]+ ".fasta"}
    df = pd.DataFrame(data=d)

    filepath = Path(output)  
    df.to_csv(filepath, index=False)

if __name__ == "__main__":
    getreference() # Runs the function