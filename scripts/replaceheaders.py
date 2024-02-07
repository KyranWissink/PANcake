def extract_fasta_headers(fasta_file):
    headers = []
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                header = line[1:].strip()  # Remove '>' and any trailing newline
                headers.append(header)
    return headers

def replace_gff_chromosomes_v1(gff_file, new_headers, output_file):
    with open(gff_file, 'r') as file:
        gff_lines = file.readlines()

    with open(output_file, 'w') as outfile:
        for line in gff_lines:
            if not line.startswith('#') and line.strip():  # Skip comments and empty lines
                parts = line.split('\t')
                scaffold_number = parts[0].split('_')[1]  # Assuming the format is 'scaffold_1', 'scaffold_2', etc.
                try:
                    scaffold_index = int(scaffold_number) - 1  # Convert to zero-based index
                    parts[0] = new_headers[scaffold_index]  # Replace chromosome name
                except (IndexError, ValueError):
                    pass  # Skip lines with unexpected format
                outfile.write('\t'.join(parts))
            else:
                outfile.write(line)

def replace_gff_chromosomes_v2(gff_file, new_headers, output_file):
    with open(gff_file, 'r') as file:
        gff_lines = file.readlines()

    header_index = 0
    chromosome_map = {}

    # Create a map for old chromosome names to new ones
    for line in gff_lines:
        if not line.startswith('#') and line.strip():
            parts = line.split('\t')
            chromosome = parts[0]
            if chromosome not in chromosome_map:
                if header_index < len(new_headers):
                    chromosome_map[chromosome] = new_headers[header_index]
                    header_index += 1
                else:
                    break  # Stop if we run out of new headers

    # Replace the chromosome names in the GFF
    with open(output_file, 'w') as outfile:
        for line in gff_lines:
            if not line.startswith('#') and line.strip():
                parts = line.split('\t')
                chromosome = parts[0]
                if chromosome in chromosome_map:
                    parts[0] = chromosome_map[chromosome]
                outfile.write('\t'.join(parts))
            else:
                outfile.write(line)

# Replace 'sequences.fa' and 'genes.gff' with your actual file paths
fasta_headers = extract_fasta_headers('sequences.fa')

# Check which method to use based on the headers
# You can add more conditions if needed
if len(fasta_headers) == 0:  # Change this condition based on your criteria
    replace_gff_chromosomes_v1('genes.gff', fasta_headers, 'updated_genes.gff')
else:
    replace_gff_chromosomes_v2('genes.gff', fasta_headers, 'updated_genes.gff')
