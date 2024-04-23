import sys

# Open the file you want to write to
with open(snakemake.log[0], 'w') as f:
    # Redirect standard output and error to the file
    sys.stdout = f
    sys.stderr = f

    import pandas as pd
    from pathlib import Path


    lineage_paths = snakemake.params[0]

    dfs = []
    for lineage, path in lineage_paths.items():
        df = pd.read_csv(path, sep='\t', low_memory=False, header=0)
        df['lineage'] = lineage
        dfs.append(df)

    df = pd.concat(dfs)
    keep_columns = ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'ID', 'Parent', 'locus', 'Name', 'description', 'old_ID', 'lineage']
    existing_columns = [column for column in keep_columns if column in df.columns]
    df = df[existing_columns]

    df.to_csv(snakemake.output[0], sep='\t', index=False)

