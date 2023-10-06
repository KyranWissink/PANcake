#################
#   FUNCTIONS   #
#################

function check_and_fill_parameters {
  # Function to check and fill missing parameters
  # Calculates or estimates parameters if missing
  # Sets default values if necessary

  if [[ -z "$number_of_genomes" ]]; then
    # Calculate the number of genomes by counting headers in the input sample
    if [[ -z "$input_sample" || ! -f "$input_sample" ]]; then
      echo "Error: Missing or invalid input sample file."
      return 1
    fi

    number_of_genomes=$(zgrep ">" "${input_sample}" | wc -l)
    echo "Missing parameter for number of genomes. Calculated: ${number_of_genomes}"
  fi

  if [[ -z "$percent_identity" ]]; then
    # Calculate percent identity using 'mash triangle'
    max_divergence=$(mash triangle -E "${input_sample}" | awk '{print $3}' | sort -g | tail -n1)
    percent_identity=$(awk "BEGIN { print 100 - $max_divergence * 100 }")
    echo "Missing parameter for percent identity. Estimated: ${percent_identity}"
    if [ "$(echo "$percent_identity < 75" | bc -l)" -eq 1 ]; then
      echo "Error: Percent identity is too low to reliably build a pangenome graph."
    fi
  fi

  if [[ -z "$poa_parameters" ]]; then
    # Determine POA parameters based on percent identity
    if [ "$(echo "$percent_identity < 99" | bc -l)" -eq 1 ]; then
      poa_parameters="asm5"
    elif [ "$(echo "$percent_identity < 90" | bc -l)" -eq 1 ]; then
      poa_parameters="asm10"
    else
      poa_parameters="asm20"
      echo "Missing parameter for poa parameters. Using ${poa_parameters} (based on percent identity)"
    fi
  fi

  if [[ -z "$segment_length" ]]; then
    segment_length=10000
    echo "Missing parameter for segment length. Using default (10000bp)"
  fi
}

function combine_fasta {
  # Function to combine multiple FASTA files into one
  # Parameters:
  #   $1: Directory containing FASTA files

  local directory="$1"
  local output_file="${directory}combined.fa"

  # Skip combination of fasta files if the combined file already exists
  if [[ -f "$output_file" ]]; then
    echo "'$output_file' already exists, skipping file combination."
    return
  fi

  # Throw an error if the directory does not exist.
  if [[ ! -d "$directory" ]]; then
    echo "Error: Directory '$directory' does not exist."
    return 1
  fi

  fasta_files=("$directory"*.fa*)

  # Throw an error if no fasta files were found in the directory.
  if [[ ${#fasta_files[@]} -eq 0 ]]; then
    echo "Error: No .fasta files found in directory '$directory'."
    return 1
  fi

  # Combine all the fasta files into one big file, and keep genome/chromosome information
  for filename in "${fasta_files[@]}"; do
    genome_name=$(basename -s .fasta "$filename" | basename -s .fa "$filename") # Cover both options
    sed "s/^>/>$genome_name\_/" "$filename" | sed -E "s/^>$genome_name>$genome_name/>$genome_name/" # Remove duplicate genome names
  done > "$output_file"

  local exit_code=$?
  if [ "$exit_code" -ne 0 ]; then
    echo "Error: Failed to combine .fasta files into '$output_file'."
    return $exit_code
  fi

  echo "Combined ${#fasta_files[@]} .fasta files into '$output_file'."
}

function run_seqpart {
  local wd="$1"

  mkdir $wd
  mv ${input_sample} ${wd}combined.fa
  bgzip -@ 4 ${wd}combined.fa
  mash dist ${wd}combined.fa.gz ${wd}combined.fa.gz -i > ${wd}distances.tsv
  python3 scripts/mash2net.py -m ${wd}distances.tsv
  python3 scripts/net2communities.py \
    -e ${wd}distances.tsv.edges.list.txt \
    -w ${wd}distances.tsv.edges.weights.txt \
    -n ${wd}distances.tsv.vertices.id2name.txt
}

function run_snakemake {
  # Function to prepare everything for snakemake to run and then run snakemake
  # Parameters:
  #   $8: Community index. If none is supplied, will assume there is no sequence partitioning.
  local number_of_genomes="$1"
  local percent_identity="$2"
  local poa_parameters="$3"
  local segment_length="$4"
  local threads="$5"
  local runid="$6"
  local input_sample="$7"
  local community="${8:-}"

  # Validate input parameters (add more validation if needed)
  if [[ -z "$number_of_genomes" || -z "$percent_identity" || -z "$poa_parameters" || -z "$segment_length" || -z "$threads" || -z "$runid" || -z "$input_sample" ]]; then
    echo "Error: Missing or empty input parameter."
    return 1
  fi

  # Define the snakemake config path using conditional assignment
  # Community number is already included in the runid.
  CONFIG_FILE="output/${runid}/snakemake_config.yaml"

  # Create a YAML-formatted config file
  cat << EOF > "$CONFIG_FILE"
pggb:
  number_of_genomes: $number_of_genomes
  percent_identity: $percent_identity
  poa_parameters: "$poa_parameters"
  segment_length: $segment_length
  threads: $threads
runid: $runid
input_sample: "$input_sample"
EOF

  # Run Snakemake with the defined configuration
  snakemake -p --forcerun --cores "$threads" --configfile "$CONFIG_FILE"

  local exit_code=$?
  if [ "$exit_code" -ne 0 ]; then
    echo "Error: Snakemake pipeline failed with exit code $exit_code."
    return $exit_code
  fi
}

function analyse_community {
  # Function to analyse a community in case of sequence partitioning
  # Parameters:
  #   $1: Community index

  local i="$1"

  if [[ -z "$i" ]]; then
    echo "Error: Missing community index."
    return 1
  fi

  echo "Analysing community: ${i}"

  local input_sample="${seqpart_dir}community.${i}.fa.gz"
  local new_runid="${runid}/community${i}"

  echo "Initialising..."
  number_of_genomes=$(zgrep ">" "${input_sample}" | wc -l)

  if [[ -z "$number_of_genomes" ]]; then
    echo "Error: Unable to determine the number of genomes."
    return 1
  fi

  echo "Running snakemake..."
  run_snakemake "$number_of_genomes" "$percent_identity" "$poa_parameters" "$segment_length" "$threads" "$new_runid" "$input_sample" "$i"

  local exit_code=$?
  if [ "$exit_code" -ne 0 ]; then
    echo "Error: Analysis for community ${i} failed with exit code $exit_code."
    return $exit_code
  fi
}
