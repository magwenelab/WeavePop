import pandas as pd
import io
import click

@click.command()
@click.option('-d', '--depth', type=click.Path(exists=True), help='Table with depth per chromosome')
@click.option('-g', '--global_mode', type=click.Path(exists=True), help='Table with global mode of sample')
@click.option('-o', '--output', type=click.Path(), help='Output file')
def normalize(depth, global_mode, output):
    depth = pd.read_csv(depth, sep='\t', header = 0)
    global_mode = pd.read_csv(global_mode, sep='\t', header = 0)
    Global_Mode = global_mode['Global_Mode'][0]
    depth['Norm_Chrom_Mean'] = depth['Chrom_Mean'] / Global_Mode
    depth['Norm_Chrom_Median'] = depth['Chrom_Median'] / Global_Mode
    depth['Norm_Global_Mean'] = depth['Global_Mean'] / Global_Mode
    depth['Norm_Global_Median'] = depth['Global_Median'] / Global_Mode
    depth = depth.round(2)
    depth.to_csv(output, sep='\t')

if __name__ == '__main__':
    normalize()