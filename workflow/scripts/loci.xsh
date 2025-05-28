#/usr/bin/env xonsh

import pandas as pd
from pathlib import Path
import click

@click.command()
@click.option('--locifile', '-l', multiple=False, required=True, type=click.Path(exists=True), help = 'Path to CSV with gene ID (locus_tag) in first column and name of locus on second column.')
@click.option('--lin_gff', '-g', required=True, type=click.Path(exists=True), help = 'Path to TSV annotation file of reference genome')
@click.option('--output', '-o', multiple=False, default = 'files/loci_to_plot.tsv', show_default=True, type=click.Path(exists=False), help = 'Path to output file.')

def getloci(locifile, lin_gff, output):
    print("Reading and concatenating annotations...")
    annotations = pd.read_csv(lin_gff, sep='\t', header=0, low_memory=False)

    print("Filtering genes...")
    level1_bool = annotations.primary_tag.str.contains("gene")     
    level1 = annotations[level1_bool]

    # if not 'locus_tag' in level1.columns:
    #     level1 = level1.assign(locus_tag = level1.ID)

    print("Reading list of gene IDs...")
    mygenes= pd.read_csv(locifile, sep=',', header=0,  names=("gene_id", "loci"))
    myloci=level1.set_index('gene_id').join(mygenes.set_index('gene_id'))
    myloci['gene_id'] = myloci.index
    myloci.dropna(subset=['loci'], inplace=True)

    print("Saving output...")
    myloci.to_csv(output,index= False, sep ='\t')

    print("Done!")

if __name__ == "__main__":
    getloci()