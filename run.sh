#!/bin/bash

set -e  # To prevent a repeat of the big scare of 04/05/2023
source functions.sh

###########################
#   INITIAL PARAMETERS    #
###########################

percent_identity=95
poa_parameters=""
segment_length=10000
threads=16
multiple_chromosomes=0

# Output colours
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[1;36m'
NC='\033[0m' # No Colour

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
    -p|--percent-identity)
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
if [[ -z "$input_sample" ]]; then
    echo -e "${GREEN}MPI Version: 0.7${NC}\n"
    echo -e "Usage:\n"
    echo -e "\tbash run.sh [options] -i <input>\n"
    echo -e "Description:\n"
    echo -e "\tPipeline for pangenome graph creation using pggb\n"
    echo -e "Options:\n"
    echo -e "Mandatory:\n"
    echo -e "\t-i --input-sample\t\tinput file(s) (fa, fa.gz, dir).\n"
    echo -e "Optional:\n"
    echo -e "\t-r --runid\t\tname for the run. Will also name directories this.\n"
    echo -e "\t-mc --multiple-chromosomes\tUse this parameter if the sample contains multiple chromosomes.\n"
    echo -e "\t-n --number-of-genomes\t\tThe number of genomes in the sample\n"
    echo -e "\t-p --percent-identity\t\tThe lowest similarity between all sequences in percentages [default: 95]\n"
    echo -e "\t-poa --poa-parameters\t\tThe partial order alignment parameters to use (asm5, asm10, asm20)\n"
    echo -e "\t-s --segment-length\t\tSegment length for mapping [default: 10k]\n"
    echo -e "\t-t --threads\t\t\tNumber of threads to use [default: 16]\n"
    exit 1
fi

# Check if optional parameters are within valid ranges or values
if [[ ! -z "$percent_identity" && ($percent_identity -lt 0 || $percent_identity -gt 100) ]]; then
  echo -e "${RED}Error${NC}: Percent identity should be between 0 and 100."
  exit 1
fi

if [[ ! -z "$poa_parameters" && ($poa_parameters != "asm5" && $poa_parameters != "asm10" && $poa_parameters != "asm20") ]]; then
  echo "${RED}Error${NC}: POA parameters should be one of asm5, asm10, or asm20."
  exit 1
fi

if [[ $segment_length -le 0 ]]; then
  echo "${RED}Error${NC}: Segment length should be above 0."
  exit 1
fi

if [[ $threads -le 0 ]]; then
  echo "${RED}Error${NC}: Threads should be above 0."
  exit 1
fi




# If no runid was supplied, then create one based on the filename
if [ -z "$runid" ]; then
  runid=$(basename "$input_sample" | cut -d. -f1)
fi

# Check if the input is a directory. If so, create a combined fasta file and adjust the sample parameter accordingly
if [[ -d "$input_sample" ]]; then
    if [[ ! $input_sample == */ ]]; then
        input_sample="${input_sample}/"
    fi
  combine_fasta $input_sample
  input_dir=${input_sample}
  input_sample=${input_dir}"combined.fa"
fi


# Call the function to check and fill missing parameters
check_and_fill_parameters

# Display the parsed parameters
echo -e "${CYAN}Sample path: ${input_sample}"
echo "Multiple Chromosomes: $multiple_chromosomes"
echo "Number of Genomes: $number_of_genomes"
echo "Percent Identity: $percent_identity"
echo "POA Parameters: $poa_parameters"
echo "Segment Length: $segment_length"
echo "Threads: $threads"
echo "RunID: $runid"
echo -e "Input Sample: ${input_sample}${NC}"



#############################
#            MAIN           #
#############################

mkdir -p "output/${runid}"

if [[ $multiple_chromosomes == 1 ]]; then

  seqpart_dir=${input_dir}seqpart/

  if [ ! -d $seqpart_dir ]; then
    run_seqpart $seqpart_dir
  else
    read -p "It looks like sequence partitioning was already performed on this data: directory ${seqpart_dir} already exists. Do you want to execute this code again? (y/n): " execute_seqpart
    if [ "$execute_seqpart" == "y" ]; then
      rm -rf $seqpart_dir
      run_seqpart $seqpart_dir
    else
      echo "Code execution skipped."
    fi
  fi


  echo "Running sequence partitioning"
  ncommunities=$(ls ${seqpart_dir} | grep distances.tsv.edges.weights.txt.community | wc -l)

  for i in $(seq 0 $(($ncommunities - 1))); do
    echo "Indexing community $i"
    mkdir -p "output/${runid}/community${i}"
    samtools faidx ${seqpart_dir}combined.fa.gz $(cat ${seqpart_dir}distances.tsv.edges.weights.txt.community.$i.txt) | \
    bgzip -@ 4 -c > ${seqpart_dir}community.$i.fa.gz
  done

  echo "Sequence partitioning finished"

  export -f analyse_community  # Export the function to make it available to parallel
  export -f run_snakemake
  export number_of_genomes
  export percent_identity
  export poa_parameters
  export segment_length
  export threads
  export runid
  export input_sample
  export input_dir
  export seqpart_dir

  for i in $(seq 0 $(($ncommunities - 1))); do
    analyse_community ${i}
  done

else  # No sequence partitioning, directly call the function
  run_snakemake $number_of_genomes $percent_identity $poa_parameters $segment_length $threads $runid $input_sample
fi

echo -e  "${GREEN}Done! Results can be found in ${PWD}/output/${runid}${NC}"


