import pandas as pd
import io
import click
from scipy import ndimage
import numpy as np

@click.command()
@click.option('-di', '--depth_input', help='Path to BED file with depth of each window.', type=click.Path(exists=True))
@click.option('-gi', '--genome_wide_depth_input', help='Path to TSV file with the genome-wide depth.', type=click.Path(exists=True))
@click.option('-do', '--depth_output', help='Path to output BED file with normalized depth.', type=click.Path())
@click.option('-s', '--smoothing_size', help='Size parameter for the smoothing function.', type=int)

def normalize(depth_input, genome_wide_depth_input, depth_output, smoothing_size):
    print("Reading depth BED file...")
    windows = pd.read_csv(depth_input, sep='\t', header=None)
    windows.columns = ['accession', 'start', 'end', 'depth']

    print("Reading genome-wide depth file...")
    genome_wide_depth_df = pd.read_csv(genome_wide_depth_input, sep='\t', header= 0)
    genome_wide_depth = genome_wide_depth_df['global_median'][0]

    print("Normalizing depth...")
    windows.loc[:,'norm_depth'] = windows['depth'] / genome_wide_depth
    windows.loc[:,'norm_depth'] = windows.loc[:,'norm_depth'].round(2)
    cov_array = np.array(windows["norm_depth"])
    smoothed_array = ndimage.median_filter(cov_array, size=smoothing_size)
    windows.loc[:,'smooth_depth']=pd.Series(smoothed_array)

    print("Saving BED file...")
    windows.to_csv(depth_output, sep='\t', index=False, header=False)

    print("Done!")

if __name__ == '__main__':
    normalize()