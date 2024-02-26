import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    
with open("logs/ploidy/test.log", "w") as f:
    sys.stderr = sys.stdout = f
    
    from scipy import ndimage
    import numpy as np
    import pandas as pd

    cov = pd.read_csv(snakemake.input[0], sep = "\t", header = 0)

    cov_array = np.array(cov["Norm_Median"])

    smoothed_array = ndimage.median_filter(cov_array, size=snakemake.params[0])

    cov['Smooth']=pd.Series(smoothed_array)

    cov = cov.loc[:, ~cov.columns.str.startswith('Global') & ~cov.columns.str.startswith('Chrom')]

    cov.to_csv(snakemake.output[0], index=False, header = True, sep='\t')