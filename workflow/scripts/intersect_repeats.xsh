import pandas as pd
import io

structure = 'results/samples/mosdepth/SRS8318899/ploidy_table.tsv'
repeats = 'results/references/VNI/repeats/05_full/VNI.bed'



intersect = $(tail -n +2 @(structure) | bedtools intersect -a stdin -b @(repeats) -wo)

df = pd.read_csv(io.StringIO(intersect), sep='\t', header=None)

header = ['sv_Accession', 'sv_Start', 'sv_End', 'sv_Size', 'sv_Norm_Cov', 'sv_Smooth_Cov', 'sv_Structure', 'r_Accession', 'r_Start', 'r_End', 'r_Type', 'Overlap_bp'] 
df.columns = header
