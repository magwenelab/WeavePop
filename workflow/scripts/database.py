import sys

log_file = snakemake.log[0]
sys.stdout = open(log_file, 'a')
sys.stderr = sys.stdout

import pandas as pd
import duckdb

metadata=snakemake.input[0]
chromosomes=snakemake.input[1]
cnvs=snakemake.input[2]
cnv_chroms=snakemake.input[3]
mapq_depth=snakemake.input[4]
gff_tsv=snakemake.input[5]
effects=snakemake.input[6]
variants=snakemake.input[7]
classif=snakemake.input[8]
presence=snakemake.input[9]
lofs=snakemake.input[10]
nmds=snakemake.input[11]
coding_sequences=snakemake.input[12]
ref_coding_sequences=snakemake.input[13]

output=snakemake.output[0]

print("Using the following arguments:")
print(f"1. metadata: {metadata}")
print(f"2. chromosomes: {chromosomes}")
print(f"3. cnvs: {cnvs}")
print(f"4. cnv_chroms: {cnv_chroms}")
print(f"5. mapq_depth: {mapq_depth}")
print(f"6. gff: {gff_tsv}")
print(f"7. effects: {effects}")
print(f"8. variants: {variants}")
print(f"9. variants classification: {classif}")
print(f"10. presence: {presence}")
print(f"11. lofs: {lofs}")
print(f"12. nmds: {nmds}")
print(f"13. coding_sequences: {coding_sequences}")
print(f"14. ref_coding_sequences: {ref_coding_sequences}")
print(f"15. output: {output}")

print("Reading metadata table...")
df_metadata = pd.read_csv(metadata) 
print("Metadata table done!")

print("Reading Copy-number variants table...")
df_cnv = pd.read_csv(cnvs, sep='\t',header = 0)
print("Copy-number variants table done!")

print("Reading Copy-number variants per chromosome table...")
df_cnv_chroms = pd.read_csv(cnv_chroms, sep='\t', dtype={'chromosome': str})
print("Copy-number variants per chromosome table done!")

print("Reading MAPQ-depth table...")
df_mapq_depth = pd.read_csv(mapq_depth, sep='\t')
print("MAPQ Depth table done!")

print("Reading chromosome table...")
df_chroms = pd.read_csv(chromosomes, header = 0, dtype={'chromosome': str, 'lenght': int})
print("Chromosome names table done!")

print("Reading GFF table...")
df_gff = pd.read_csv(gff_tsv, sep='\t', header = 0, low_memory=False)
print("GFF table done!")

print("Reading effects table...")
df_effects = pd.read_csv(effects, header = 0, sep='\t')
print("Effects table done!")

print("Reading variants table...")
df_variants = pd.read_csv(variants, header = 0, sep='\t')
print("Variants table done!")

print("Reading variant classification table...")
df_classif = pd.read_csv(classif, header = 0, sep='\t')
print("Variant classification table done!")

print("Reading presence table...")
df_presence = pd.read_csv(presence, header = 0, sep='\t')
print("Presence table done!")

print("Reading lofs table...")
df_lofs = pd.read_csv(lofs, header = 0, sep='\t')
print("Loss of function table done!")

print("Reading nonsense-mediated decay table...")
df_nmds = pd.read_csv(nmds, header = 0, sep='\t')
print("Nonsense-mediated decay table done!")

print("Reading coding_sequences table...")
df_coding_sequences = pd.read_csv(coding_sequences, header = 0, sep=',')
print("Coding sequences table done!")

print("Reading reference coding coding_sequences table...")
df_ref_coding_sequences = pd.read_csv(ref_coding_sequences, header = 0, sep=',')
print("Reference coding_sequences table done!")

print("Formatting dataframes done!")

print("Create function to add table to db")
def create_table_with_constraints(con, df, table_name, constraints):
    """
    Create a DuckDB table from a DataFrame with specified constraints.

    Args:
        con: duckdb.DuckDBPyConnection object
        df: pandas.DataFrame to insert
        table_name: str, name of the final table
        constraints: list of str, each a constraint (e.g. "PRIMARY KEY (id)", "FOREIGN KEY (user_id) REFERENCES users(user_id)")
    """
    # Register DataFrame as a DuckDB view
    view_name = f"df_{table_name}"
    con.register(view_name, df)

    # Get column definitions from DataFrame dtypes
    dtype_map = {
        'int64': 'BIGINT',
        'float64': 'DOUBLE',
        'object': 'VARCHAR',
        'bool': 'BOOLEAN'
    }
    columns = []
    for col, dtype in df.dtypes.items():
        duck_type = dtype_map.get(str(dtype), 'VARCHAR')
        columns.append(f'"{col}" {duck_type}')
    # Add constraints
    columns += constraints
    create_stmt = f"CREATE TABLE {table_name} ({', '.join(columns)})"
    con.execute(create_stmt)
    # Insert data
    con.execute(f"INSERT INTO {table_name} SELECT * FROM {view_name}")

print("Connecting to database")
con = duckdb.connect(database=output)

print("Adding dataframes to database")

create_table_with_constraints(con, df_chroms, 'chromosomes',[
    "PRIMARY KEY (accession)",
])

create_table_with_constraints(con, df_metadata, 'metadata', [
    "PRIMARY KEY (sample)"
])
create_table_with_constraints(con, df_cnv, 'cnvs',[
    'PRIMARY KEY (accession,start,"end",sample)',
    "FOREIGN KEY (accession) REFERENCES chromosomes(accession)",
    "FOREIGN KEY (sample) REFERENCES metadata(sample)"
])
create_table_with_constraints(con, df_cnv_chroms, 'cnv_chroms',[
    "PRIMARY KEY (sample,accession,cnv)",
    "FOREIGN KEY (accession) REFERENCES chromosomes(accession)",
    "FOREIGN KEY (sample) REFERENCES metadata(sample)"
])
create_table_with_constraints(con, df_gff, 'gff', [
    "PRIMARY KEY (feature_id,accession)",
    "FOREIGN KEY (accession) REFERENCES chromosomes(accession)",
])
create_table_with_constraints(con, df_mapq_depth, 'mapq_depth',[
    "PRIMARY KEY (sample,feature_id)",
    "FOREIGN KEY (sample) REFERENCES metadata(sample)",
])
create_table_with_constraints(con, df_coding_sequences, 'coding_sequences', [
    "PRIMARY KEY (sample,feature_id,seq_type)",
    "FOREIGN KEY (sample) REFERENCES metadata(sample)"
])
create_table_with_constraints(con, df_ref_coding_sequences, 'ref_coding_sequences', [
    "PRIMARY KEY (ref_genome,feature_id,seq_type)"
])
create_table_with_constraints(con, df_variants, 'variants', [
    "PRIMARY KEY (var_id)",
    "FOREIGN KEY (accession) REFERENCES chromosomes(accession)"
])
create_table_with_constraints(con, df_classif, 'variant_classification', [
    "PRIMARY KEY (var_id)",
    "FOREIGN KEY (var_id) REFERENCES variants(var_id)"
])
create_table_with_constraints(con, df_presence, 'presence',[
    "PRIMARY KEY (var_id,sample)",
    "FOREIGN KEY (var_id) REFERENCES variants(var_id)",
    "FOREIGN KEY (sample) REFERENCES metadata(sample)"
])
create_table_with_constraints(con, df_lofs, 'lofs',[
    "PRIMARY KEY (var_id,gene_name)",
    "FOREIGN KEY (var_id) REFERENCES variants(var_id)"
])
create_table_with_constraints(con, df_nmds, 'nmds',[
    "PRIMARY KEY (var_id,gene_name)",
    "FOREIGN KEY (var_id) REFERENCES variants(var_id)"
])
create_table_with_constraints(con, df_effects, 'effects',[
    "FOREIGN KEY (var_id) REFERENCES variants(var_id)"
])

print("Closing connection to database")
con.close()
print("Done!")