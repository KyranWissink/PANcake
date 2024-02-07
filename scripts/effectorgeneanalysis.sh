#!/bin/bash

set -e

# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 ID_EXTRACT_SCRIPT NODE_TO_GENES_SCRIPT"
  exit 1
fi

# Define paths to your data and scripts
ID_EXTRACT_SCRIPT=$1
NODE_TO_GENES_SCRIPT=$2
INDEX_GIRAFFE="index.giraffe.gbz"
INDEX_XG="index.xg"

# Iterate over each community folder
for community_dir in ./community*; do
    echo "Processing $community_dir..."
    cd $community_dir
    # Perform vg annotate
    vg annotate -x index.xg -f ~/workflow/software/snpEff/data/combined.gff -t 64 > annotated.gam

    # Convert gam to json
    vg view -aj "annotated.gam" > "annotated.json"

    # Extract node ids
    python3 $ID_EXTRACT_SCRIPT < "annotated.json" > "annotated_ids.txt"

    # Apply the process to each Folac directory
    for folac_dir in Folac*; do
        echo "Processing $folac_dir..."
        cd $folac_dir
        vg pack -x ../index.xg -i "aln.pack" -N "../annotated_ids.txt" -d -Q 20 -t 90 | awk '$4 > 10' | awk '{print $2}' | sort -u > "annotated_nodes.txt"

        # Convert the nodes to genes
        python3 $NODE_TO_GENES_SCRIPT
        cd ..
    done
    cd ..
done

