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
if [[ -z "$runid" || -z "$input_sample" ]]; then
    echo -e "${GREEN}MPI Version: 0.7${NC}\n"
    echo -e "Usage:\n"
    echo -e "\tbash run.sh [options] -i <input> -r <runid>\n"
    echo -e "Description:\n"
    echo -e "\tPipeline for pangenome graph creation using pggb\n"
    echo -e "Options:\n"
    echo -e "Mandatory:\n"
    echo -e "\t-i --input-sample\t\tinput fasta file or directory. Zipped files are supported.\n"
    echo -e "\t-r --runid\t\tdirectory to output all files. Need not be yet created.\n"
    echo -e "Optional:\n"
    echo -e "\t-mc --multiple-chromosomes\tUse this parameter if the sample contains multiple chromosomes.\n"
    echo -e "\t-n --number-of-genomes\t\tThe number of genomes in the sample\n"
    echo -e "\t-p --percent-identity\t\tThe lowest similarity between all sequences in percentages [default: 95]\n"
    echo -e "\t-poa --poa-parameters\t\tThe partial order alignment parameters to use (asm5, asm10, asm20)\n"
    echo -e "\t-s --segment-length\t\tSegment length for mapping [default: 10k]\n"
    echo -e "\t-t --threads\t\t\tNumber of threads to use [default: 16]\n"
    exit 1
fi


# Check if the input is a directory. If so, create a combined fasta file and adjust the sample parameter accordingly
echo "Sample path: ${input_sample}"
runid=$(basename "$input_sample")

if [[ -d "$input_sample" ]]; then
    if [[ ! input_sample == */ ]]; then
        input_sample="${input_sample}/"
    fi
  combine_fasta $input_sample
  input_dir=${input_sample}
  input_sample=${input_sample}"combined.fa"
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

# Call the function to check and fill missing parameters
check_and_fill_parameters

# Display the parsed parameters
echo -e "${CYAN}Multiple Chromosomes: $multiple_chromosomes"
echo "Number of Genomes: $number_of_genomes"
echo "Percent Identity: $percent_identity"
echo "POA Parameters: $poa_parameters"
echo "Segment Length: $segment_length"
echo "Threads: $threads"
echo "RunID: $runid"
echo "Input Sample: ${input_sample}${NC}"

#############################
#            MAIN           #
#############################

if [[ $multiple_chromosomes == 1 ]]; then
  echo "Running sequence partitioning"

  if [ ! -d ${input_dir}seqpart/ ]; then
    mkdir ${input_dir}seqpart/
    mv ${input_sample} ${input_dir}seqpart/combined.fa
    wd=${input_dir}seqpart/
    bgzip -@ 4 ${wd}combined.fa

    mash dist ${wd}combined.fa.gz ${wd}combined.fa.gz -i > ${wd}distances.tsv
    python3 scripts/mash2net.py -m ${wd}distances.tsv
    python3 scripts/net2communities.py \
      -e ${wd}distances.tsv.edges.list.txt \
      -w ${wd}distances.tsv.edges.weights.txt \
      -n ${wd}distances.tsv.vertices.id2name.txt
  else
    wd=${input_dir}seqpart/
  fi

  echo "Indexing data"
  ncommunities=$(ls ${wd} | grep distances.tsv.edges.weights.txt.community | wc -l)

  for i in $(seq 0 $(($ncommunities - 1))); do
    echo "Indexing community $i"
    samtools faidx ${wd}combined.fa.gz $(cat ${wd}distances.tsv.edges.weights.txt.community.$i.txt) | \
    bgzip -@ 4 -c > ${wd}community.$i.fa.gz
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
  export wd

  for i in $(seq 0 $(($ncommunities - 1))); do
    echo "Analysing community ${i}"
    analyse_community ${i}
  done

else  # No sequence partitioning, directly call the function
  run_snakemake $number_of_genomes $percent_identity $poa_parameters $segment_length $threads $runid $input_sample
fi

echo "${GREEN}Done! Results can be found in ${PWD}/output/${runid}${NC}"


