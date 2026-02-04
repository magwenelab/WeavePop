import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import pandas as pd

lin_gff=snakemake.input[0]
output=snakemake.output[0]


if len(snakemake.input) > 1:
    locifile=snakemake.input[1]

    print("Reading and concatenating annotations...")
    annotations = pd.read_csv(lin_gff, sep='\t', header=0, low_memory=False)

    print("Filtering genes...")
    level1_bool = annotations.primary_tag.str.contains("gene")     
    level1 = annotations[level1_bool]

    print("Reading list of gene IDs...")
    mygenes= pd.read_csv(locifile, sep=',', header=0,  names=("gene_id", "loci"))
    myloci=level1.set_index('gene_id').join(mygenes.set_index('gene_id'))
    myloci['gene_id'] = myloci.index
    myloci.dropna(subset=['loci'], inplace=True)

    print("Saving output...")
    myloci.to_csv(output,index= False, sep ='\t')

    print("Done!")

else:
    myloci = pd.DataFrame(columns=['gene_id', 'loci'])
    myloci.to_csv(output,index= False, sep ='\t')
    print("No loci file provided, created empty output.")
    print("Done!")