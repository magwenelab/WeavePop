#/usr/bin/env xonsh

import pandas as pd
from pathlib import Path
import click

@click.command()
@click.option('--genefile', '-g', multiple=False, required=True, type=click.Path(exists=True), help = 'Path to CSV with gene ID (locus_tag) in first column and name of locus on second column.')
@click.option('--output', '-o', multiple=False, default = 'files/loci_to_plot.tsv', show_default=True, type=click.Path(exists=False), help = 'Path to output file.')
@click.argument('referencetsv', nargs=-1, required=True, type=click.Path(exists=True)) # Path to TSV annotation file of reference genome

def getloci(genefile, referencetsv, output):
    """This script creates an annotation table <output> with only level 1 features of the provided gene IDs (locus_tag) in <genefile>.
    A column 'loci' with the corresponding locus name is added.
    The output table has the annotation of the genes in all the reference genomes REFERENCETSV given as positional arguments."""
    print("Reading and concatenating annotations...")
    mydfs = []
    for lin in referencetsv:
        mydfs.append(pd.read_csv(Path(lin), sep='\t', header=0, low_memory=False))
    annotations = pd.concat(mydfs)

    print("Filtering genes...")
    level1_bool = annotations.primary_tag.str.contains("gene")     
    level1 = annotations[level1_bool]

    # if not 'locus_tag' in level1.columns:
    #     level1 = level1.assign(locus_tag = level1.ID)

    print("Reading list of gene IDs...")
    mygenes= pd.read_csv(Path(genefile), sep=',', header=0,  names=("gene_id", "loci"))
    myloci=level1.set_index('gene_id').join(mygenes.set_index('gene_id'))
    myloci['gene_id'] = myloci.index
    myloci.dropna(subset=['loci'], inplace=True)

    print("Saving output...")
    myloci.to_csv(output,index= False, sep ='\t')

    print("Done!")

if __name__ == "__main__":
    getloci()