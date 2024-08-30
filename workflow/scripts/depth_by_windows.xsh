import pandas as pd
import io
import click
from scipy import ndimage
import numpy as np

@click.command()
@click.option('-di', '--depth_input', help='Path to BED file with depth of each window.', type=click.Path(exists=True))
@click.option('-gi', '--global_mode_input', help='Path to TSV file woth the mode of the global depth.', type=click.Path(exists=True))
@click.option('-do', '--depth_output', help='Path to output BED file with normalized depth.', type=click.Path())
@click.option('-s', '--smoothing_size', help='Size parameter for the smoothing function.', type=int)

def normalize(depth_input, global_mode_input, depth_output, smoothing_size):
    print("Read depth BED file.")
    windows = pd.read_csv(depth_input, sep='\t', header=None)
    windows.columns = ['Accession', 'Start', 'End', 'Depth']

    print("Read global mode file.")
    global_mode = pd.read_csv(global_mode_input, sep='\t', header= 0)
    Global_Mode = global_mode['Global_Mode'][0]

    print("Normalize depth.")
    windows.loc[:,'Norm_Depth'] = windows['Depth'] / Global_Mode
    windows.loc[:,'Norm_Depth'] = windows.loc[:,'Norm_Depth'].round(2)
    cov_array = np.array(windows["Norm_Depth"])
    smoothed_array = ndimage.median_filter(cov_array, size=smoothing_size)
    windows.loc[:,'Smooth_Depth']=pd.Series(smoothed_array)
    windows.to_csv(depth_output, sep='\t', index=False, header=False)

if __name__ == '__main__':
    normalize()