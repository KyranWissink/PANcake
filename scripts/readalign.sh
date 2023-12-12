#!/bin/bash

set -e

# Default values
output_dir="readmapping"
threads=64

# Parse command-line options
while getopts ":o:t:f:i:" opt; do
  case $opt in
    o) output_dir="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    f)
      if [ -z "$fastq_file1" ]; then
        fastq_file1="$OPTARG"
      elif [ -z "$fastq_file2" ]; then
        fastq_file2="$OPTARG"
      else
        echo "Error: Too many -f options. Please provide only two files for paired-end reads." >&2
        exit 1
      fi
      ;;
    i) input_dir="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

# Check if input directory is provided
if [ -z "$input_dir" ]; then
  echo "Error: Input directory (-i) is required."
  exit 1
fi

# Check if at least one input fastq file is provided
if [ -z "$fastq_file1" ] && [ -z "$fastq_file2" ]; then
  echo "Error: At least one input fastq file is required."
  exit 1
fi

# Iterate through communities in the input directory
for community_dir in "$input_dir"/*/; do
  community=$(basename "$community_dir")
  community_output_dir="$community/$output_dir"

  mkdir -p "$community_output_dir"  # Create readmapping directory for each community

  # Perform your operations for each community
  echo "Indexing $community"
  cd $community_output_dir
  vg autoindex -w giraffe -g "../data.gfa" -t "$threads" -R XG
  echo "Aligning reads"

  if [ -n "$fastq_file1" ] && [ -n "$fastq_file2" ]; then
    # Paired-end reads
    echo "Using paired-end reads"
    vg giraffe -x index.xg -Z index.giraffe.gbz -m index.min -d index.dist -f "$fastq_file1" -f "$fastq_file2" -t "$threads" -o GAM > "aln.gam"
  elif [ -n "$fastq_file1" ]; then
    # Single-end reads
    echo "Using single-end reads"
    vg giraffe -x index.xg -Z index.giraffe.gbz -m index.min -d index.dist -f "$fastq_file1" -t "$threads" -o GAM > "aln.gam"
  fi

  echo "Packing"
  vg pack -x index.xg -g "aln.gam" -Q 5 -o "aln.pack"

  echo "Calling"
  vg call index.xg -k "aln.pack" -a -d 1 -t "$threads" > "vg_calls.vcf"

  cd ../..  # Return to the main directory
  sleep 5
done
