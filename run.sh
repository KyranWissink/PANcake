#!/bin/bash

set -e  # To prevent a repeat of the big scare of 04/05/2023

###########################
#   INITIAL PARAMETERS    #
###########################

percent_identity=""
poa_parameters=""
segment_length=10000
threads=16
multiple_chromosomes=0

#################
#   FUNCTIONS   #
#################

combine_fasta() {
    local directory="$1"
    local output_file="${directory}combined.fa"

    if [[ ! $directory == */ ]]; then
        directory="$directory/"
    fi

    if [[ -f $output_file ]]; then
        echo "'$output_file' already exists, skipping file combination."
        return
    fi

    if [[ ! -d $directory ]]; then
        echo "Directory '$directory' does not exist."
        exit 1
    fi

    fasta_files=(${directory}*.fa*)

    if [[ ${#fasta_files[@]} -eq 0 ]]; then
        echo "No .fasta files found in directory '$directory'."
        exit 1
    fi

    for filename in "${fasta_files[@]}"; do
        file_name=$(basename "$filename" | cut -d. -f1)
        sed "s/^>/>$file_name_/" "$filename"
    done > "$output_file"

    echo "Combined ${#fasta_files[@]} .fasta files into '$output_file'."
}

# Function to run necessary processes
function run_snakemake {

  # Define the path for the YAML config file
  CONFIG_FILE="snakemake_config.yaml"

  # Create a YAML-formatted config file
  echo "multiple_chromosomes: $multiple_chromosomes" > "$CONFIG_FILE"
  echo "pggb:" >> "$CONFIG_FILE"
  echo "  number_of_genomes: $number_of_genomes" >> "$CONFIG_FILE"
  echo "  percent_identity: $percent_identity" >> "$CONFIG_FILE"
  echo "  poa_parameters: $poa_parameters" >> "$CONFIG_FILE"
  echo "  segment_length: $segment_length" >> "$CONFIG_FILE"
  echo "  threads: $threads" >> "$CONFIG_FILE"
  echo "runid: $runid" >> "$CONFIG_FILE"
  echo "input_sample: $input_sample" >> "$CONFIG_FILE"

  # ... (rest of your script)

  # Run Snakemake with the defined configuration
  snakemake -p --forcerun --cores --configfile "$CONFIG_FILE"

}

# Function to analyze a single community
function analyze_community {
  local i=$1
  echo "Analysing community: ${i}"
  sed "s#sample:.*#sample: ${wd}community.${i}.fa.gz#g" config.yaml > temp.yaml && mv temp.yaml config.yaml
  sed "s#runid:.*#runid: ${runid}/community${i}#g" config.yaml > temp.yaml && mv temp.yaml config.yaml

  # Call the function to run processes
  run_snakemake
}



#########################
#   PARAMETER PARSING   #
#########################

# Parse command line options
while [[ $# -gt 0 ]]; do
  case "$1" in
    -mc|--multiple-chromosomes)
      multiple_chromosomes=1
      shift
      ;;
    -n|--number-of-genomes)
      number_of_genomes="$2"
      shift 2
      ;;
    -pi|--percent-identity)
      percent_identity="$2"
      shift 2
      ;;
    -poa|--poa-parameters)
      poa_parameters="$2"
      shift 2
      ;;
    -s|--segment-length)
      segment_length="$2"
      shift 2
      ;;
    -t|--threads)
      threads="$2"
      shift 2
      ;;
    -r|--runid)
      runid="$2"
      shift 2
      ;;
    -i|--input-sample)
      input_sample="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Check for required parameters
if [[ -z "$runid" || -z "$input_sample" ]]; then
    echo -e "Version: 0.7\n"
    echo -e "Usage:\n"
    echo -e "\tbash run.sh [options] -i <input> -o <output_dir>\n"
    echo -e "Description:\n"
    echo -e "\tPipeline for pangenome graph creation using pggb\n"
    echo -e "Options:\n"
    echo -e "Mandatory:\n"
    echo -e "\t-i --input-sample\t\tinput fasta file or directory. Zipped files are supported.\n"
    echo -e "\t-o --output-directory\t\tdirectory to output all files. Need not be yet created.\n"
    echo -e "Optional:\n"
    echo -e "\t-mc --multiple-chromosomes\tUse this parameter if the sample contains multiple chromosomes.\n"
    echo -e "\t-n --number-of-genomes\t\tThe number of genomes in the sample\n"
    echo -e "\t-p --percent-identity\t\tThe lowest similarity between all sequences in percentages\n"
    echo -e "\t-poa --poa-parameters\t\tThe partial order alignment parameters to use (asm5, asm10, asm20)\n"
    echo -e "\t-s --segment-length\t\tSegment length for mapping [default: 10k]\n"
    echo -e "\t-t --threads\t\t\tNumber of threads to use [default: 16]\n"
    exit 1
fi


# Check if the input is a directory. If so, create a combined fasta file and adjust the sample parameter accordingly
[[ -d "$sample" ]] && sample=$(combine_fasta "$sample")

# Check if optional parameters are within valid ranges or values
if [[ ! -z "$percent_identity" && ($percent_identity -lt 0 || $percent_identity -gt 100) ]]; then
  echo "Error: Percent identity should be between 0 and 100."
  exit 1
fi

if [[ ! -z "$poa_parameters" && ($poa_parameters != "asm5" && $poa_parameters != "asm10" && $poa_parameters != "asm20") ]]; then
  echo "Error: POA parameters should be one of asm5, asm10, or asm20."
  exit 1
fi

if [[ $segment_length -le 0 ]]; then
  echo "Error: Segment length should be above 0."
  exit 1
fi

if [[ $threads -le 0 ]]; then
  echo "Error: Threads should be above 0."
  exit 1
fi

# Check missing optional parameters
if [[ -z "$number_of_genomes" ]]; then
  # Calculate the number of genomes by counting headers in the input sample
  number_of_genomes=$(zgrep ">" "${input_sample}" | wc -l)
  echo "Missing parameter for number of genomes. Calculated: ${number_of_genomes}"
fi

if [[ -z "$percent_identity" ]]; then
  # Calculate percent identity using 'mash triangle'
  max_divergence=$(mash triangle -E "${input_sample}" | awk '{print $3}' | sort -g | tail -n1)
  percent_identity=$(awk "BEGIN { print 100 - $max_divergence * 100 }")
  echo "Missing parameter for percent identity. Estimated: ${lowest_percent_identity}"
fi

if [[ -z "$poa_parameters" ]]; then
  # Determine POA parameters based on percent identity
  if (( percent_identity > 99 )); then
    poa_parameters="asm5"
  elif (( percent_identity > 90 )); then
    poa_parameters="asm10"
  else
    poa_parameters="asm20"
    echo "Missing parameter for poa parameters. Using ${poa_parameters} (based on percent identity)"
  fi
fi


if [[ -z "$segment_length" ]]; then
  segment_length=10000
  echo "Missing parameter for segment length. Using default (3000bp)"
fi


# Display the parsed parameters
echo "Multiple Chromosomes: $multiple_chromosomes"
echo "Number of Genomes: $number_of_genomes"
echo "Percent Identity: $percent_identity"
echo "POA Parameters: $poa_parameters"
echo "Segment Length: $segment_length"
echo "Threads: $threads"
echo "RunID: $runid"
echo "Input Sample: $input_sample"

#############################
#            MAIN           #
#############################

if [[ $multiple_chromosomes == 1 ]]; then
  echo "Running sequence partitioning"

  if [ ! -d ${sample}seqpart/ ]; then
    python3 scripts/combine.py ${sample}
    mkdir ${sample}seqpart/
    mv ${sample}combined.fa ${sample}/seqpart/combined.fa
    wd=${sample}seqpart/
    bgzip -@ 4 ${wd}combined.fa

    mash dist ${wd}combined.fa.gz ${wd}combined.fa.gz -i > ${wd}distances.tsv
    python3 scripts/mash2net.py -m ${wd}distances.tsv
    python3 scripts/net2communities.py \
      -e ${wd}distances.tsv.edges.list.txt \
      -w ${wd}distances.tsv.edges.weights.txt \
      -n ${wd}distances.tsv.vertices.id2name.txt
  else
    wd=${sample}seqpart/
  fi

  echo "Indexing data"
  ncommunities=$(ls ${wd} | grep distances.tsv.edges.weights.txt.community | wc -l)

  # Use parallel to analyze communities in parallel
  parallel -j $threads 'analyze_community {}' ::: $(seq 0 $ncommunities)

  echo "Sequence partitioning finished."

else  # No sequence partitioning, directly call the function
  run_snakemake
fi

echo "Done! Results can be found in ${PWD}/output/${runid}"


