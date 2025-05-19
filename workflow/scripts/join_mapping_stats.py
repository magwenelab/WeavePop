import sys

log_file=snakemake.log[0]
sys.stdout = open(log_file, "a")
sys.stderr = sys.stdout

import pandas as pd

input=snakemake.input
output=snakemake.output[0]
min_depth=snakemake.params.min_depth
min_high_mapq=snakemake.params.min_high_mapq
min_pm=snakemake.params.min_pm
min_coverage=snakemake.params.min_coverage

print("Reading and concatenating files...")
stats = pd.concat([pd.read_csv(f, sep="\t") for f in input])

print(f"Minimum depth percentage: {min_depth}")
print(f"Minimum high MAPQ percentage: {min_high_mapq}")
print(f"Minimum mapped percentage: {min_pm}")
print(f"Minimum coverage percentage: {min_coverage}")

print("Adding quality warning flag...")
stats['mapq_warning'] = stats.apply(lambda row: "MAPQ-Low" if row['percent_high_mapq'] < min_high_mapq else None, axis=1)
stats['pm_warning'] = stats.apply(lambda row: "Mapped-Low" if row['percent_mapped'] < min_pm else None, axis=1)
stats['depth_warning'] = stats.apply(lambda row: "Depth-Low" if row['genome_median_depth_good'] < min_depth else None, axis=1)
stats['coverage_warning'] = stats.apply(lambda row: "Coverage-Low" if row['coverage_good'] < min_coverage else None, axis=1)

print("Joining warnings...")
stats['quality_warning'] = stats[['mapq_warning', 'pm_warning', 'depth_warning','coverage_warning']].apply(lambda x: '_'.join(x.dropna()), axis=1)
stats = stats.drop(columns = ['mapq_warning', 'pm_warning', 'depth_warning','coverage_warning'])

print("Quality warnings added:")
print(stats[['sample','quality_warning']])

print("Saving mapping stats table...")
stats.to_csv(output, index=False, sep = "\t")
print("Done!")
