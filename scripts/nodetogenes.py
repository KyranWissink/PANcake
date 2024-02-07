import csv
import json
import sys

# Increase the maximum field size limit
csv.field_size_limit(sys.maxsize)

# Load node lengths from the CSV file
node_lengths = {}
with open('../nodes.csv', 'r') as csv_file:
    reader = csv.reader(csv_file, delimiter=';')
    next(reader)  # Skip header
    for row in reader:
        node_id, sequence = row
        node_lengths[node_id] = len(sequence)

# Load the node IDs
with open('nodes.txt') as f:
    node_ids = {line.strip() for line in f}

# Parse the JSON file and calculate aligned nucleotide count
gene_to_info = {}
with open('../annotated.json', 'r') as file:
    for line in file:
        alignment = json.loads(line)
        gene_id = alignment.get('name', 'Unknown')

        # Initialize data structure for each gene
        if gene_id not in gene_to_info:
            gene_to_info[gene_id] = {'aligned_length': 0, 'node_count': set()}

        # Check each mapping in the path
        if 'path' in alignment and 'mapping' in alignment['path']:
            for mapping in alignment['path']['mapping']:
                if 'position' in mapping and 'node_id' in mapping['position']:
                    node_id = str(mapping['position']['node_id'])  # Convert to string for matching

                    if node_id in node_ids:
                        total_edit_length = sum(edit.get('to_length', 0) for edit in mapping['edit'])
                        gene_to_info[gene_id]['aligned_length'] += total_edit_length
                        gene_to_info[gene_id]['node_count'].add(node_id)

# Output the mapping in TSV format
with open('gene_id_node_count_aligned_length.tsv', 'w') as out_file:
    out_file.write('GeneID\tNodeCount\tAlignedLength\n')
    for gene_id, info in gene_to_info.items():
        out_file.write(f'{gene_id}\t{len(info["node_count"])}\t{info["aligned_length"]}\n')

print("TSV file created: gene_id_node_count_aligned_length.tsv")
