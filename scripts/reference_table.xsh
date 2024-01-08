from collections import defaultdict
from pathlib import Path
import sys
import pandas as pd
import click


# Create dataframe from SAMPLESTABLE and REFGENOMETABLE with the sample, read filenames, lineage and reference assembly file information
@click.command()
@click.option('-s', '--samplestable', 'samplestable', required=True, type=str, help = "Sample metadata table with columns sample and group")
@click.option('-r', '--refgenometable', 'refgenometable',required=True, type=str, help = "Table with reference groups (lineages) and their corresponding fasta files, with columns group and file, respectively")
@click.option('-o', '--output', 'output', required=True, type=str, help = "File name of output table with the samples and the corresponding reference genome files.")
@click.option('-f1', '--suffix1', 'suffix1', required=True, type=str, help = "Suffix of the forward fastq filenames.")
@click.option('-f2', '--suffix2', 'suffix2', required=True, type=str, help = "Suffix of the referse fastq filenames.")

def getreference(samplestable,refgenometable,output, suffix1, suffix2):
    lineage = pd.read_csv(samplestable)
    lineage = lineage[["sample","group"]]

    d={'sample': lineage["sample"], 'file1': lineage["sample"]+ suffix1, 'file2': lineage["sample"]+ suffix2}
    df = pd.DataFrame(data=d)
    df = df.set_index('sample').join(lineage.set_index('sample'))
    df.reset_index(inplace=True)

    reference = pd.read_csv(refgenometable)
    reference = reference[["group","file"]]
    reference = reference.rename(columns={'file': 'refgenome'})

    df = df.set_index('group').join(reference.set_index('group'))
    df.reset_index(inplace=True)

    filepath = Path(output)  
    df.to_csv(filepath, index=False)

if __name__ == "__main__":
    getreference() # Runs the function