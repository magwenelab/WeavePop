import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import pandas as pd
import numpy as np
from scipy import ndimage

depth_input=snakemake.input.depth
depth_output=snakemake.output.windows
smoothing_size=snakemake.params.smoothing_size
sample=snakemake.wildcards.sample

print("Reading depth BED file...")
windows = pd.read_csv(depth_input, sep='\t', header=None)
windows.columns = ['accession', 'start', 'end', 'depth']

print("Calculating genome-wide median depth from windows...")
genome_median_depth = windows['depth'].median().round(4)

print("Normalizing depth...")
windows.loc[:,'norm_depth'] = windows['depth'] / genome_median_depth
windows.loc[:,'norm_depth'] = windows.loc[:,'norm_depth'].round(2)
cov_array = np.array(windows["norm_depth"])
smoothed_array = ndimage.median_filter(cov_array, size=smoothing_size)
windows.loc[:,'smooth_depth']=pd.Series(smoothed_array)

print("Saving BED file...")
windows.to_csv(depth_output, sep='\t', index=False, header=False)

print("Done!")