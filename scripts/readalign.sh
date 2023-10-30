#!/bin/bash

# Define the list of communities
communities=("community0" "community1" "community2" "community3" "community4" "community5" "community6" "community7" "community8" "community9" "community10" "community11" "community12")

# Loop through each community
for community in "${communities[@]}"; do
    mkdir -p "$community/readmapping"  # Create readmapping directory for each community
    cd "$community/readmapping" || exit 1

    # Perform your operations for each community
    echo "Indexing"
    vg autoindex -w giraffe -g "../data.gfa" -t 64
    echo "Aligning reads"
    vg giraffe -Z index.giraffe.gbz -m index.min -d index.dist -f "../../TR4_II5_R1_merged.fastq" -f "../../TR4_II5_R2_merged.fastq" -t 64 -o bam > aln.bam
    echo "Sorting"
    samtools sort -o sorted.bam aln.bam
    echo "Calculating coverage"
    samtools depth -a sorted.bam > coverage.txt
    awk -F'\t' '{sum[$1] += $3; count[$1]++} END {OFS="\t"; for (path in sum) print "Path:", path, "Average:", sum[path]/count[path]}' coverage.txt > average_coverage.txt

    cd ../..  # Return to the main directory
done
