import click
import pandas as pd

@click.command()
@click.option('--input','-i', help='Input file: TSV version of a GFF file', required=True, type = click.Path(exists=True))
@click.option('--output_tsv', '-ot', help='Output TSV', required=True, type = click.Path())
@click.option('--output_gff', '-og', help='Output GFF', required=True, type = click.Path())

def fix_gff(input, output_tsv, output_gff):
    df = pd.read_csv(input, sep='\t', header=[0])
    keep_columns = ['seq_id', 'source_tag', 'primary_tag', 'start', 'end', 'score', 'strand', 'frame', 'ID', 'Parent', 'locus', 'Name', 'description']
    exisiting_columns = [column for column in keep_columns if column in df.columns]
    df = df[exisiting_columns]
    df['source_tag'] = 'DiversityPipeline'
    df['primary_tag'] = df['primary_tag'].replace('protein_coding_gene', 'gene')
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
    df_fixed['old_ID'] = df_fixed['ID']
    df_fixed['ID'] = df_fixed['new_ID']
    df_fixed.drop(['feature_number', 'suffix', 'new_ID'], axis=1, inplace=True)

    print("Saving TSV")
    df_fixed.to_csv(output_tsv, sep='\t', index=False, header=True, na_rep='')

    print("Making GFF")
    attribute_columns = ['ID', 'Parent', 'locus', 'Name', 'description', 'old_ID', 'old_Parent']
    existing_attributes = [column for column in attribute_columns if column in df_fixed.columns]
    df_gff = df_fixed.copy()
    for column in existing_attributes:
        df_gff[column] = df_gff[column].apply(lambda x: column + '=' + x if pd.notnull(x) else x)
    df_gff['attributes'] = df_gff[existing_attributes].apply(lambda x: ';'.join(x.dropna()), axis=1)
    df_gff.drop(existing_attributes, axis=1, inplace=True)
    df_gff['strand'] = df_gff['strand'].replace(-1,'-').replace(1,'+')

    print("Saving GFF")
    df_gff.to_csv(output_gff, sep='\t', index=False, header=False)

    
if __name__ == '__main__':
    fix_gff()