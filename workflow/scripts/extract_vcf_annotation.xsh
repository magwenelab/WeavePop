import pandas as pd
import click
import io
import os
import re
import vcf

@click.command()
@click.option("--vcf_path", "-i", type = click.Path(exists=True), help="Path to the VCF file")
@click.option("--gff_tsv", "-g", type = click.Path(), help="Path to the table with annotation of the references")
@click.option("--effects_out", "-e", type = click.Path(), help="Path to the effects output file")
@click.option("--variants_out", "-v", type = click.Path(), help="Path to the variants output file")
@click.option("--lofs_out", "-f", type = click.Path(), help="Path to the lofs output file")
@click.option("--nmds_out", "-n", type = click.Path(), help="Path to the nmds output file")
@click.option("--lineage", "-l", type = str, help="Lineage")

def extract_annotation(vcf_path, gff_tsv, effects_out, variants_out, lofs_out, nmds_out, lineage):
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
            var_id = var_id[0]
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
                transcript_id = effect_info_parts[8]
                exon_rank = effect_info_parts[9]
                data_effects.append([var_id,  effect_type, impact, effect, codon_change, aa_change, aa_length, gene_name, biotype, coding, transcript_id, exon_rank])
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
    df_variants = df_variants[['var_id', 'lineage', 'accession', 'pos', 'ref', 'alt']]

    print("LOFs dataframe")
    df_lofs = pd.DataFrame(data_lofs, columns=['var_id', 'gene_name', 'num_transcripts', 'percent_affected'])
    print("NMDs dataframe")
    df_nmds = pd.DataFrame(data_nmds, columns=['var_id', 'gene_name', 'num_transcripts', 'percent_affected'])

    print("Effects dataframe")
    effects_pre = pd.DataFrame(data_effects, columns=['var_id', 'effect_type', 'impact','effect', 'codon_change', 'amino_acid_change', 'amino_acid_length', 'gene_name', 'transcript_biotype', 'gene_coding', 'transcript_id', 'exon_rank'])

    print("Formating effects dataframe and adding gene IDs from reference GFF...")
    print("Reading GFF file...")
    df_gff = pd.read_csv(gff_tsv, sep='\t', header = 0, low_memory=False)

    print("Getting gene IDs from GFF file...")  
    id_cols = ['locus', 'Name', 'ID']
    existing_id_cols = [col for col in id_cols if col in df_gff.columns]
    gff_ids = df_gff[existing_id_cols].drop_duplicates().copy()
    rename_dict = {'locus': 'gene_id', 'Name': 'gene_name', 'ID': 'feature_id'}
    existing_columns = gff_ids.columns
    filtered_rename_dict = {k: v for k, v in rename_dict.items() if k in existing_columns}
    gff_ids.rename(columns=filtered_rename_dict, inplace=True)
    gff_ids = gff_ids.dropna(subset=['gene_id'])

    print("Subsetting effects table...")

    effects_pre.replace('', pd.NA, inplace=True)
    df_gene_transcript = effects_pre[(effects_pre['gene_name'].notnull()) & (effects_pre['transcript_id'].notnull())].copy()
    df_gene_no_transcript = effects_pre[(effects_pre['gene_name'].notnull()) & (effects_pre['transcript_id'].isnull())].copy()
    df_no_gene_no_transcript = effects_pre[(effects_pre['gene_name'].isnull()) & (effects_pre['transcript_id'].isnull())].copy()

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
    df_gene_transcript['gene_id'] = df_gene_transcript['transcript_id'].map(feature_to_gene)
    print("Concatenating dataframes...")
    df_effects = pd.concat([df_gene_transcript, df_gene_no_transcript_fixed, df_no_gene_no_transcript])
    print("Finished formatting effects dataframe!")

    print("Saving dataframes to CSV files...")
    df_variants.to_csv(variants_out, sep = "\t",  index=False)
    df_effects.to_csv(effects_out, sep = "\t", index=False)
    df_lofs.to_csv(lofs_out, sep = "\t", index=False)
    df_nmds.to_csv(nmds_out, sep = "\t",index=False)

    print("Done!")

if __name__ == '__main__':
    extract_annotation()