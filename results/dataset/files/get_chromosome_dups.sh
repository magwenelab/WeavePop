#!/bin/bash

input_file="structural_variants.tsv"
output_file="filtered_variants.tsv"

# Extracting the header line from the input file
header=$(head -n 1 "$input_file")

# Filtering out four NC1 large chr1 duplications, plus odd sample with strange plotting
grep "VNI_CP003820.1" "$input_file" | grep "Duplication" | grep -vE "PMY3622|PMY3624|PMY3626|PMY3618|PMY3457" | grep -v "Repetitive sequence" | sort -k2 -n > "$output_file"

# Prepending the header line to the output file
echo "$header" | cat - "$output_file" > temp && mv temp "$output_file"

cut -f 1 "$input_file" | sort -u > chromosomes.tsv