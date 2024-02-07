#!/bin/bash

set -e

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 ALIGN_DIR OUTPUT"
  exit 1
fi

# Directory where trimmed alignment files are stored
ALIGN_DIR=$1
SUPERMATRIX=$2

# Initialize an associative array to hold the concatenated sequences for each species
declare -A supermatrix_sequences

# Loop through each alignment file
for file in "$ALIGN_DIR"/*.fa; do
    # Extract sequences from the file and append them to the corresponding species in the associative array
    while read -r line; do
        if [[ $line == ">"* ]]; then
            current_species=$line
        else
            supermatrix_sequences[$current_species]+=$line
        fi
    done < "$file"
done

# Write the concatenated sequences to the supermatrix file
for species in "${!supermatrix_sequences[@]}"; do
    echo "$species" >> "$SUPERMATRIX"
    echo "${supermatrix_sequences[$species]}" >> "$SUPERMATRIX"
done

