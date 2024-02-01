#!/usr/bin/env python


"""

Parse a FASTA file, writing (appending) each record to a file named with that 
record's ID.  Each record is prepended with the given sample name .

For example, the FASTA file:

```
>gene01
ATGC
>gene02
TTGCA
```

when called with  the sample name "Sample1" will produce:

1. `gene01.fasta`

```
>Sample1_gene01
ATGC
```

2. `gene02.fasta`

```
>Sample1_gene02
ATGC
```    

"""

import os, fcntl
from pathlib import Path

import click
from Bio import SeqIO




@click.command()
@click.option("--outdir",
              type=click.Path(dir_okay=True, file_okay=False),
              default="./",
              help="Directory to write output FASTA files to.")
@click.argument("input", type=click.Path(dir_okay=False,exists=True))
@click.argument("sample", type=str)
def output_by_id(input, sample, outdir):
    """
    Parse a FASTA file (INPUT), writing  each record to a file named with that 
    record's ID.  The sample name (SAMPLE) is added to the header of each record.

    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for record in SeqIO.parse(input, "fasta"):
        old_id = record.id 
        new_id = f"{sample}_{old_id}"
        record.description = f"{sample}_{record.description}"
        record.id = new_id
        with open(Path(outdir) / f"{old_id}.fasta", "a") as output_handle:
            fcntl.flock(output_handle, fcntl.LOCK_EX)
            SeqIO.write(record, output_handle, "fasta")
            fcntl.flock(output_handle, fcntl.LOCK_UN)


if __name__ == '__main__':
    output_by_id()