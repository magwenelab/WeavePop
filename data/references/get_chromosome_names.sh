#!/bin/bash

output_file="/FastData/kwilhoit/DiversityPipeline/data/references/headers.csv"
fasta_dir="/FastData/kwilhoit/DiversityPipeline/data/references"

# Remove existing output file if it exists
if [ -f "$output_file" ]; then
    rm "$output_file"
fi

# Loop through each .fasta file in the directory
for fasta_file in "$fasta_dir"/*.fasta; do
    # Get the filename without the extension
    filename=$(basename "$fasta_file" .fasta)
    
    # Initialize the counter
    counter=1
    
    # Read each line of the .fasta file
    while IFS= read -r line; do
        # Check if the line starts with a '>'
        if [[ $line == ">"* ]]; then
            # Extract the header without the starting '>'
            header=$(echo "$line" | cut -c2-)
            
            # Extract the first part of the header before the space
            column1=$(echo "$header" | cut -d' ' -f1)
            
            # Extract the first part of the string in column 1 before the underscore
            column2=$(echo "$column1" | cut -d'_' -f1)
            
            # Write the columns to the output file with column 1 and 2 switched
            echo "$column2,$column1,$counter" >> "$output_file"
            
            # Increment the counter
            counter=$((counter + 1))
            
            # Reset the counter if it reaches 15
            if [ $counter -eq 15 ]; then
                counter=1
            fi
        fi
    done < "$fasta_file"
done