#!/bin/bash

# Set the input and output file paths
input_file="vni_samples_only_metadata.csv"
output_file="subsample_metadata.csv"

# Get the first row from the input file
head -n 1 "$input_file" > "$output_file"

# Get a random subset of 5 rows from the input file (excluding the first row)
tail -n +2 "$input_file" | shuf -n 5 >> "$output_file"