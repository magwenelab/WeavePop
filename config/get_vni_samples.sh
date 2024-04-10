#!/bin/bash

input_file="/FastData/kwilhoit/DiversityPipeline/config/no_hybrids_sample_metadata.csv"
output_file="/FastData/kwilhoit/DiversityPipeline/config/vni_samples_only_metadata.csv"

# Create an empty output file
> "$output_file"
# Copy the first row from the input file to the output file
head -n 1 "$input_file" >> "$output_file"
# Read each line from the input file
while IFS= read -r line; do
    # Extract the first column value
    first_column=$(echo "$line" | cut -d ',' -f 1)

    # Check if the first column value matches "VNI" exactly
    if [ "$first_column" = "VNI" ]; then
        # Append the matching line to the output file
        echo "$line" >> "$output_file"
    fi
done < "$input_file"