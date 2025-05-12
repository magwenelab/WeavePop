import pandas as pd
import io
import click
from scipy import ndimage
import numpy as np

@click.command()
@click.option('-di', '--depth_input', help='Path to BED file with depth of each window.', type=click.Path(exists=True))
@click.option('-do', '--depth_output', help='Path to output BED file with normalized depth.', type=click.Path())
@click.option('-ss', '--smoothing_size', help='Size parameter for the smoothing function.', type=int)
@click.option('-sm', '--sample', help='Sample name.', type=str)


def normalize(depth_input, depth_output, smoothing_size, sample):
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

if __name__ == '__main__':
    normalize()