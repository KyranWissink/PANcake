#!/bin/bash

# Set the parent directory where genome folders are located
PARENT_DIR="$1"

# Check if the parent directory is provided
if [ -z "$PARENT_DIR" ]; then
    echo "Usage: $0 PARENT_DIR"
    exit 1
fi

# Create an associative array to hold the counts of each ortholog
declare -A ortholog_counts

# Initialize a variable to count the number of genomes
genome_count=0

# Loop through each genome directory
for genome_dir in "$PARENT_DIR"/*busco_output/run_ascomycota_odb10/busco_sequences/single_copy_busco_sequences; do
    # Check if the genome directory exists and is a directory
    if [ -d "$genome_dir" ]; then
        # Increment genome count
        ((genome_count++))

        # Find all ortholog .fa files and loop through them
        for ortholog in "$genome_dir"/*.faa; do
            # Check if the ortholog file exists and is a file
            if [ -f "$ortholog" ]; then
                # Extract the name of the ortholog
                ortholog_name=$(basename "$ortholog" .fa)

                # Increment the count for this ortholog
                if [ -n "${ortholog_counts[$ortholog_name]}" ]; then
                    ortholog_counts[$ortholog_name]=$(( ${ortholog_counts[$ortholog_name]} + 1 ))
                else
                    ortholog_counts[$ortholog_name]=1
                fi
            fi
        done
    fi
done

# Now, filter orthologs present in all genomes
output_file="$PARENT_DIR/orthologs.txt"
for ortholog in "${!ortholog_counts[@]}"; do
    if [[ ${ortholog_counts[$ortholog]} -eq $genome_count ]]; then
        echo "$ortholog" >> "$output_file"
    fi
done

echo "Orthologs present in all genomes are listed in $output_file"
