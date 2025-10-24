import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import pandas as pd
import vcf

vcf_path = snakemake.input.vcf
gff_tsv = snakemake.input.tsv
metadata = snakemake.input.metadata
presence = snakemake.input.presence
variants_out = snakemake.output.variants
effects_out = snakemake.output.effects
lofs_out = snakemake.output.lofs
nmds_out = snakemake.output.nmds
classif_out = snakemake.output.classif
lineage = snakemake.wildcards.lineage

print("Getting tables from SnpEff result...")
data_effects = []
data_variants = []
data_lofs = []
data_nmds = []
print("Iterating over VCF file records to extract information...")
with open(vcf_path, 'r') as file:
    # Open the VCF file
    vcf_reader = vcf.Reader(file)
    # Add a unique identifier to each record
    for i, record in enumerate(vcf_reader):
        var_id = record.INFO['var_id']
        data_variants.append([var_id, record.CHROM, record.POS, record.REF, record.ALT])
        # Get the SnpEff annotations from the INFO field
        annotations = record.INFO.get('EFF', [])
        # Get the LOF annotations from the INFO field
        lofs = record.INFO.get('LOF', [])
        # Get the NMD annotations from the INFO field
        nmds = record.INFO.get('NMD', [])
        # Iterate over all annotations
        for annotation in annotations:
            # Split the annotation into its components
            annotation_parts = annotation.split('(')
            effect_type = annotation_parts[0]
            effect_info = annotation_parts[1]
            effect_info_parts = effect_info.split('|')
            impact = effect_info_parts[0]
            effect = effect_info_parts[1]
            codon_change = effect_info_parts[2]
            aa_change = effect_info_parts[3]
            aa_length = effect_info_parts[4]
            gene_name = effect_info_parts[5]
            biotype = effect_info_parts[6]
            coding = effect_info_parts[7]
            feature_id = effect_info_parts[8]
            exon_rank = effect_info_parts[9]
            data_effects.append([var_id,  effect_type, impact, effect, codon_change, aa_change, aa_length, gene_name, biotype, coding, feature_id, exon_rank])
        # Iterate over all lofs
        for lof in lofs:
            lof_parts = lof.split('|')
            first = lof_parts[0]
            gene_name = first.replace('(', '')
            num_transcripts = lof_parts[2]
            last = lof_parts[3]
            percent_affected = last.replace(')', '')
            data_lofs.append([var_id, gene_name, num_transcripts, percent_affected])
        # Iterate over all nmds
        for nmd in nmds:
            nmd_parts = nmd.split('|')
            first = nmd_parts[0]
            gene_name = first.replace('(', '')
            num_transcripts = nmd_parts[2]
            last = nmd_parts[3]
            percent_affected = last.replace(')', '')
            data_nmds.append([var_id, gene_name, num_transcripts, percent_affected])

print("Creating dataframes...")
print("Variants dataframe")
df_variants = pd.DataFrame(data_variants, columns=['var_id', 'accession', 'pos', 'ref', 'alt'])
df_variants['alt'] = df_variants['alt'].astype(str)
df_variants['alt'] = df_variants['alt'].str.replace('[', '').str.replace(']', '')
df_variants['lineage'] = lineage
df_variants = df_variants[['var_id', 'accession', 'pos', 'ref', 'alt', 'lineage']]

print("LOFs dataframe")
df_lofs = pd.DataFrame(data_lofs, columns=['var_id', 'gene_name', 'num_transcripts', 'percent_affected'])
print("NMDs dataframe")
df_nmds = pd.DataFrame(data_nmds, columns=['var_id', 'gene_name', 'num_transcripts', 'percent_affected'])

print("Effects dataframe")
effects_pre = pd.DataFrame(data_effects, columns=['var_id', 'effect_type', 'impact','effect', 'codon_change', 'amino_acid_change', 'amino_acid_length', 'gene_name', 'transcript_biotype', 'gene_coding', 'feature_id', 'exon_rank'])

print("Formating effects dataframe and adding gene IDs from reference GFF...")
print("Reading GFF file...")
df_gff = pd.read_csv(gff_tsv, sep='\t', header = 0, low_memory=False)

print("Getting gene IDs from GFF file...")  
gff_ids = df_gff[['gene_id', 'gene_name', 'feature_id']].drop_duplicates().copy()
gff_ids = gff_ids.dropna(subset=['gene_id'])

print("Subsetting effects table...")

effects_pre.replace('', pd.NA, inplace=True)
df_gene_transcript = effects_pre[(effects_pre['gene_name'].notnull()) & (effects_pre['feature_id'].notnull())].copy()
df_gene_no_transcript = effects_pre[(effects_pre['gene_name'].notnull()) & (effects_pre['feature_id'].isnull())].copy()
df_no_gene_no_transcript = effects_pre[(effects_pre['gene_name'].isnull()) & (effects_pre['feature_id'].isnull())].copy()

print("Defining function to replace gene names with gene IDs...")
def replace_with_gene_id(part):
    if part in gff_ids['gene_id'].values:
        new_part = part
    elif part in gff_ids['gene_name'].values:
        new_part = gff_ids.loc[gff_ids['gene_name'] == part, 'gene_id'].values[0]
    return new_part

print("Getting variant effects with fused gene names...")
df_fused_genes = df_gene_no_transcript[df_gene_no_transcript['gene_name'].str.contains('\\+')].copy()
if df_fused_genes.shape[0] > 0:
    print("Fused gene names found")
    print("Separating fused gene names...")
    df_fused_genes[['part1', 'part2']] = df_fused_genes['gene_name'].str.split('+', expand=True)
    print("Replacing gene names with gene IDs in part1...")
    df_fused_genes['gene_tag_id1'] = df_fused_genes['part1'].apply(replace_with_gene_id)
    print("Replacing gene names with gene IDs in part2...")
    df_fused_genes.loc[df_fused_genes['part2'].notnull(), 'gene_tag_id2'] = df_fused_genes.loc[df_fused_genes['part2'].notnull(), 'part2'].apply(replace_with_gene_id)
    print("Joining part1 with part2...")
    df_fused_genes['gene_id'] = df_fused_genes.apply(lambda row: row['gene_tag_id1'] + '+' + row['gene_tag_id2'] if pd.notnull(row['part2']) else row['gene_tag_id1'], axis=1)
    print("Removing unnecessary columns...")
    df_fused_genes.drop(['part1', 'part2', 'gene_tag_id1', 'gene_tag_id2'], axis=1, inplace=True)
    df_gene_no_transcript_fixed = df_fused_genes.copy()
else:
    print("No fused gene names found")
    df_gene_no_transcript_fixed = df_gene_no_transcript.copy()
    df_gene_no_transcript_fixed['gene_id'] = df_gene_no_transcript_fixed['gene_name'].apply(replace_with_gene_id)

print("Getting unique gene IDs from GFF...")
gff_ids_unique = gff_ids.drop_duplicates(subset='feature_id', keep='first')
print("Creating dictionary to map feature IDs to gene IDs...")
feature_to_gene = gff_ids_unique.set_index('feature_id')['gene_id']
print("Mapping transcript IDs to gene IDs...")
df_gene_transcript['gene_id'] = df_gene_transcript['feature_id'].map(feature_to_gene)
print("Concatenating dataframes...")
df_effects = pd.concat([df_gene_transcript, df_gene_no_transcript_fixed, df_no_gene_no_transcript])
print("Finished formatting effects dataframe!")

print("Classifying variants according to privateness")
df_presence = pd.read_csv(presence, sep='\t', header = 0, low_memory=False)
df_private = df_presence.groupby(['var_id']).size().reset_index(name = "n_samples")

df_metadata = pd.read_csv(metadata, header = 0)
df_samples_in_lineage= df_metadata.groupby(['lineage']).size().reset_index(name='samples_in_lineage')
df_private['samples_in_lineage'] = df_samples_in_lineage.loc[df_samples_in_lineage['lineage'] == lineage, ['samples_in_lineage']]['samples_in_lineage'].iloc[0]
  
df_private['category'] = df_private.apply(lambda row: "Reference private" if row['n_samples'] == row['samples_in_lineage']
                                                            else "Private" if row['n_samples'] == 1
                                                            else "Non-private", axis = 1)
df_private = df_private[['var_id', 'category']]


print("Classifying variants according to impact")

df_impact = df_effects.groupby(['var_id','impact']).size().reset_index(name='n_effects')
df_impact['one_effect'] = 1

df_impact = df_impact.pivot(index='var_id', columns='impact' , values='one_effect').reset_index()
df_impact = df_impact.fillna(0)

df_impact['impact'] = df_impact.apply(lambda row: "High" if row['HIGH'] == 1 
                                      else "Moderate" if row['MODERATE'] == 1 
                                      else "Low" if row['LOW'] == 1 
                                      else "Modifier", axis =1)

df_impact = df_impact[['var_id', 'impact']]

print("Merging both classifications into one table")
df_classif = pd.merge(df_private, df_impact, on = 'var_id', how = 'outer')

print("Classification table")
print(df_classif)

print("Saving dataframes to CSV files...")
df_variants.to_csv(variants_out, sep = "\t",  index=False)
df_effects.to_csv(effects_out, sep = "\t", index=False)
df_lofs.to_csv(lofs_out, sep = "\t", index=False)
df_nmds.to_csv(nmds_out, sep = "\t",index=False)
df_classif.to_csv(classif_out, sep = "\t",index=False)

print("Done!")