#/usr/bin/env xonsh

import pandas as pd
from pathlib import Path
import click

@click.command()
@click.option('--genefile', '-g', multiple=True, required=True, type=click.Path(exists=True), help = 'Path to text file with list of gene IDs (each in new line) of the locus of interest, file name will be used as locus name.')
@click.option('--output', '-o', multiple=False, default = 'files/loci_interest.tsv', show_default=True, type=click.Path(exists=False), help = 'Path to output file.')
@click.argument('referencetsv', nargs=-1, required=True, type=click.Path(exists=True)) # Path to TSV annotation file of reference genome

def getloci(genefile, referencetsv, output):
    """This script creates an annotation table <output> with the columns: 
    "seq_id", "primary_tag", "start", "description", "gene_id", "ID", "Name", "Loci" for each gene in <genefile>. 
    It takes as many gene files as desired, each with the IDs of the genes in a locus of interest.
    The output table has the annotation of the genes in all the reference genomes REFERENCETSV given as positional arguments."""
    mydfs = []
    for lin in referencetsv:
        mydfs.append(pd.read_csv(Path(lin), sep='\t', header=0, low_memory=False))
    annotations = pd.concat(mydfs)
    level1_bool = annotations.primary_tag.str.contains("gene")     
    level1 = annotations[level1_bool]

    mylocus = []
    for locus in genefile:
        name = locus.split('/')[-1].split('.')[-2]
        gene_list = list(pd.read_csv(Path(locus), header = None)[0])
        gene = level1[level1.apply(lambda row: any(value in row.values for value in gene_list), axis=1)]
        gene = gene.assign(Loci=name)
        mylocus.append(gene)
    all_loci = pd.concat(mylocus)

    all_loci.to_csv(output,index= False, sep ='\t')

if __name__ == "__main__":
    getloci()