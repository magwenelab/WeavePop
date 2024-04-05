# agat_convert_sp_gxf2gxf.pl -g data/main_reference/FungiDB-65_CneoformansH99.gff -o H99_temp.gff
# agat_sq_add_locus_tag.pl --gff H99_temp.gff --li ID -o H99_temp2.gff
# agat_sp_manage_attributes.pl --gff H99_temp2.gff --tag product/description -o H99.gff
# agat_convert_sp_gff2tsv.pl --gff H99.gff -o H99.tsv

# agat_convert_sp_gxf2gxf.pl -g data/references/North.gff -o North_temp.gff
# agat_sq_add_locus_tag.pl --gff North_temp.gff --li ID -o North_temp2.gff
# agat_sp_manage_attributes.pl --gff North_temp2.gff --tag product/description -o North.gff
# agat_convert_sp_gff2tsv.pl --gff North.gff -o North.tsv

import click
import pandas as pd

@click.command()
@click.option('--input','-i', help='Input file: TSV version of a GFF file', required=True, type = click.Path(exists=True))
@click.option('--output_tsv', '-ot', help='Output TSV', required=True, type = click.Path())
@click.option('--output_gff', '-og', help='Output GFF', required=True, type = click.Path())

def fix_gff(input, output_tsv, output_gff):
    df = pd.read_csv(input, sep='\t', header=[0])
    keep_columns = ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'ID', 'Parent', 'locus', 'Name', 'description']
    df = df[keep_columns]
    df['source_tag'] = 'DiversityPipeline'
    print("Defining feature number and suffix")
    df['feature_number'] = df.groupby(['primary_tag', 'Parent']).cumcount()+ 1
    df['feature_number'] = df['feature_number'].astype(float).fillna(0).astype(int)
    df['suffix'] = df['primary_tag'].str[:2]

    print("Defining new IDs and Parents")
    for index, row in df.iterrows():
        if row['ID'] == row['locus']:
            df.at[index, 'new_ID'] = row['ID']
        elif row['Parent'] == row['locus']:
            df.at[index, 'new_ID'] = row['locus'] + '-' + row['suffix'] + '_' + str(row['feature_number'])
        else:       
            new_Parent = df.loc[df['ID'] == row['Parent'], 'new_ID'].values[0]
            new_ID = new_Parent + '-' + row['suffix'] + str(row['feature_number'])
            df.at[index, 'new_ID'] = new_ID
            df.at[index, 'Parent'] = new_Parent
            
    print("Fixing columns")
    df_fixed = df.copy()
    df_fixed['ID'] = df_fixed['new_ID']
    df_fixed.drop(['feature_number', 'suffix', 'new_ID'], axis=1, inplace=True)
    
    print("Saving TSV")
    df_tsv = df_fixed.copy()
    df_tsv.drop(['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame'], axis=1, inplace=True)
    df_tsv.to_csv(output_tsv, sep='\t', index=False, header=True)

    print("Saving GFF")
    df_gff = df_fixed.copy()
    df_gff.drop(['Parent', 'locus', 'Name', 'description'], axis=1, inplace=True)
    df_gff['ID'] = 'ID=' + df_gff['ID'].astype(str)
    df_gff.to_csv(output_gff, sep='\t', index=False, header=False)

if __name__ == '__main__':
    fix_gff()