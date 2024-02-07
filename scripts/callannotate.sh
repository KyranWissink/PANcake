#!/bin/bash

set -e

# Check if three arguments are provided
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 VCF_FILE_PATTERN ANNOTATION_REFERENCE SNPEFF_JAR_PATH"
  exit 1
fi

# Assign arguments to variables
vcf_pattern=$1
annotation_ref=$2
SNPEFF_JAR_PATH=$3

# Default value for threads
threads=64

# Iterate through community directories
for community_dir in community*; do
  community=$(basename "$community_dir")
  echo "Processing community: $community_dir"

  # Iterate through Folac directories within the community
  for folac_dir in "$community_dir"/Folac*; do
    folac=$(basename "$folac_dir")
    echo "Processing Folac: $folac_dir"

    # Change to Folac directory
    cd "$folac_dir"

    # Get all variables
    ref_path=$(grep "$vcf_pattern" ../data.gfa | awk '{print $2}' | head -n1)
    call_out="${ref_path}.vcf"
    annotated_out="${ref_path}_annotated.vcf"
    trimmed_out="${ref_path}_annotated_trimmed.vcf"



    # If the ref_path is non-empty
    if [[ -n "$ref_path" ]]; then
      if [ ! -s "$call_out" ]  && [ ! -s "$trimmed_out" ] && [ ! -s "$annotated_out" ]; then
        echo "Calling $ref_path"
        vg call index.xg -k aln.pack -a -d 1 -t "$threads" -p "$ref_path" > "$call_out"
      fi


      # If annotated file does not exist
      if  [ ! -s "$annotated_out" ]  && [ ! -s "$trimmed_out" ]; then
        echo "Annotating $call_out"
        java -jar -Xmx64G "$SNPEFF_JAR_PATH" ann "$annotation_ref" "$call_out" > "$annotated_out"
      fi

      # If trimmed file does not exist
      if [ ! -s "$trimmed_out" ]; then
        echo "Trimming $annotated_out"
        # Trim for intergenic regions and save as trimmed file
        grep -v "intergenic_region" "$annotated_out" > "$trimmed_out"

      fi

      # Delete original and annotated VCF files
      if [ -f "$call_out" ]; then
        rm -f "$call_out"
      fi

      if [ -f "$annotated_out" ]; then
        rm -f "$annotated_out"
      fi

      echo "Processed $ref_path in $folac_dir"
    # Change back to community directory
    else
      echo "Skipping $community; no path found"
    fi
    cd ..
  done

  # Change back to the output directory
  cd ..
done

echo "Process completed."
