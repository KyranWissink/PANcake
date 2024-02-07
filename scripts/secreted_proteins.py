from Bio import SeqIO
import sys

# Check if three arguments are provided
if len(sys.argv) != 3:
    print("Usage: {} signalp_output_file proteins_file".format(sys.argv[0]))
    sys.exit(1)

# Assigning command line arguments to variables
signalp_output_file = sys.argv[1]
proteins_file = sys.argv[2]

# Store secreted protein IDs
secreted_proteins = set()

# Read secreted protein IDs from SignalP output file
with open(signalp_output_file, 'r') as file:
    for line in file:
        # Skip header lines and empty lines
        if line.startswith('#') or not line.strip():
            continue

        # Split the line into columns
        parts = line.split()
        protein_id = parts[0]
        prediction = parts[1]

        # Check if the protein is predicted to be secreted (SP(Sec/SPI))
        if prediction == 'SP(Sec/SPI)':
            secreted_proteins.add(protein_id)

# Minimum length for protein sequences
MIN_LENGTH = 10

# Filter protein sequences and write secreted proteins to a new file
with open(proteins_file, 'r') as proteins, open('secreted_proteins.fasta', 'w') as secreted_file:
    for record in SeqIO.parse(proteins, 'fasta'):
        if record.id in secreted_proteins and len(record.seq) >= MIN_LENGTH:
            SeqIO.write(record, secreted_file, 'fasta')
