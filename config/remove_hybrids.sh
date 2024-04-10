#!/bin/bash

input_file="Montoya_Corrected_sample_metadata.csv"
output_file="no_hybrids.csv"

# Loop through each row in the input file
while IFS= read -r line; do
    # Check if the first column is empty or contains an underscore
    if [[ -z ${line%%,*} || $line == *_* ]]; then
        continue  # Skip the row if it's empty or contains an underscore
    fi

    # Append the row to the output file
    echo "$line" >> "$output_file"
done < "$input_file"