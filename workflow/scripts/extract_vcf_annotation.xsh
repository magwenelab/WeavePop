import pandas as pd
import click
import io
import os
import re
import vcf

@click.command()
@click.option("--vcf_path", "-i", type = click.Path(exists=True), help="Path to the VCF file")
@click.option("--effects_out", "-e", type = click.Path(), help="Path to the effects output file")
@click.option("--variants_out", "-v", type = click.Path(), help="Path to the variants output file")
@click.option("--lofs_out", "-f", type = click.Path(), help="Path to the lofs output file")
@click.option("--nmds_out", "-n", type = click.Path(), help="Path to the nmds output file")
@click.option("--lineage", "-l", type = str, help="Lineage")

def extract_annotation(vcf_path, effects_out, variants_out, lofs_out, nmds_out, lineage):
    print("Getting tables from SnpEff result")
    data_effects = []
    data_variants = []
    data_lofs = []
    data_nmds = []
    # Iterate over the records in the VCF file
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
                locus = effect_info_parts[5]
                biotype = effect_info_parts[6]
                coding = effect_info_parts[7]
                transcript_id = effect_info_parts[8]
                exon_rank = effect_info_parts[9]
                data_effects.append([var_id,  effect_type, impact, effect, codon_change, aa_change, aa_length, locus, biotype, coding, transcript_id, exon_rank])
            # Iterate over all lofs
            for lof in lofs:
                lof_parts = lof.split('|')
                first = lof_parts[0]
                locus = first.replace('(', '')
                nb_transcripts = lof_parts[2]
                last = lof_parts[3]
                percent_transcripts = last.replace(')', '')
                data_lofs.append([var_id, locus, nb_transcripts, percent_transcripts])
            # Iterate over all nmds
            for nmd in nmds:
                nmd_parts = nmd.split('|')
                first = nmd_parts[0]
                locus = first.replace('(', '')
                nb_transcripts = nmd_parts[2]
                last = nmd_parts[3]
                percent_transcripts = last.replace(')', '')
                data_nmds.append([var_id, locus, nb_transcripts, percent_transcripts])

    print("Creating dataframes")
    df_variants = pd.DataFrame(data_variants, columns=['var_id', 'accession', 'pos', 'ref', 'alt'])
    df_effects = pd.DataFrame(data_effects, columns=['var_id', 'effect_type', 'impact','effect', 'codon_change', 'amino_acid_change', 'amino_acid_length', 'locus', 'transcript_biotype', 'gene_coding', 'transcript_id', 'exon_rank'])
    df_effects['effect_id'] = 'eff_' + lineage + '_' + (df_effects.index + 1).astype(str)
    df_effects.insert(0, 'effect_id', df_effects.pop('effect_id'))
    df_lofs = pd.DataFrame(data_lofs, columns=['var_id', 'locus', 'nb_transcripts', 'percent_transcripts'])
    df_nmds = pd.DataFrame(data_nmds, columns=['var_id', 'locus', 'nb_transcripts', 'percent_transcripts'])

    print("Adding lineage to dataframes")
    df_variants['lineage'] = lineage
    df_effects['lineage'] = lineage
    df_lofs['lineage'] = lineage
    df_nmds['lineage'] = lineage

    print("Saving dataframes to CSV files")
    df_variants.to_csv(variants_out, index=False)
    df_effects.to_csv(effects_out, index=False)
    df_lofs.to_csv(lofs_out, index=False)
    df_nmds.to_csv(nmds_out, index=False)

if __name__ == '__main__':
    extract_annotation()