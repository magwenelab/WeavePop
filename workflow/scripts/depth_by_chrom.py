import pandas as pd
import io
import click

@click.command()
@click.option('-b', '--depth_bed', help='Path to BED file with depth of each region.', type=click.Path(exists=True))
@click.option('-g', '--global_mode', help='Path to file with genome-wide mode of depth.', type=click.Path(exists=True))
@click.option('-o', '--output', help='Path to output table of stats by chromosome.', type=click.Path())
@click.option('-s', '--sample_name', help='Sample name as a string.', type=str)
def depth_by_chrom(depth_bed, output, sample_name, global_mode=None):
    print("Importing data")
    region_depth = pd.read_csv(depth_bed, sep='\t')
    region_depth.columns = ['Accession', 'Start', 'End', 'Depth']
    print("Calculating chromosome mean and median")
    Chrom_Mean = region_depth.groupby('Accession').agg(Chrom_Mean=('Depth', 'mean')).reset_index()
    Chrom_Median = region_depth.groupby('Accession').agg(Chrom_Median=('Depth', 'median')).reset_index()
    Global_Mean = region_depth['Depth'].mean()
    Global_Median = region_depth['Depth'].median()
    chrom_stats = pd.merge(Chrom_Mean, Chrom_Median, on='Accession')
    chrom_stats['Global_Mean'] = Global_Mean
    chrom_stats['Global_Median'] = Global_Median
    chrom_stats['Sample'] = sample_name
    chrom_stats = chrom_stats.round(2)
    # Put Sample column before the rest
    chrom_stats = chrom_stats[['Sample', 'Accession', 'Chrom_Mean', 'Chrom_Median', 'Global_Mean', 'Global_Median']]
    
    if global_mode:
        print("Importing global mode")
        df_global_mode = pd.read_csv(global_mode, sep='\t', header = 0)
        print("Adding global mode")
        chrom_stats['Global_Mode'] = df_global_mode['Global_Mode'][0]
        print("Normalizing chromosome depth")
        chrom_stats['Norm_Chrom_Mean'] = chrom_stats['Chrom_Mean'] / chrom_stats['Global_Mode']
        chrom_stats['Norm_Chrom_Median'] = chrom_stats['Chrom_Median'] / chrom_stats['Global_Mode']
        chrom_stats['Norm_Global_Mean'] = chrom_stats['Global_Mean'] / chrom_stats['Global_Mode']
        chrom_stats['Norm_Global_Median'] = chrom_stats['Global_Median'] / chrom_stats['Global_Mode']
        chrom_stats = chrom_stats.round(2)
        
    print("Writing output")
    chrom_stats.to_csv(output, sep='\t', index=False)

if __name__ == '__main__':
    depth_by_chrom()