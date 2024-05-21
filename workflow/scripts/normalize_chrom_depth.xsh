import pandas as pd
import io
import click

@click.command()
@click.option('-d', '--depth', type=click.Path(exists=True), help='Table with depth per chromosome')
@click.option('-g', '--global_mode', type=click.Path(exists=True), help='Table with global mode of sample')
@click.option('-s', '--sample', type=str, help='Sample name')
@click.option('-o', '--output', type=click.Path(), help='Output file')
def normalize(depth, global_mode, sample, output):


if __name__ == '__main__':
    normalize()