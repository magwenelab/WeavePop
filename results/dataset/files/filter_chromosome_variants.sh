#!/bin/bash

input_file="structural_variants.tsv"
chromosomes_file="chromosomes.tsv"

# Get the first line of the input file
first_line=$(head -n 1 "$input_file")

# Loop over each chromosome value in chromosomes.tsv
while IFS= read -r chromosome; do
    output_file="structural_variants_by_chromosome/${chromosome}_variants.tsv"  # Generate the output file name
    
    # Append the first line to the output file
    echo "$first_line" > "$output_file"
    
    # Filter and sort the data, excluding the first line
    grep "$chromosome" "$input_file" | grep "Duplication" | grep -vE "PMY3457" | grep -v "Repetitive sequence" | sort -k13 > "$output_file"
done < "$chromosomes_file"