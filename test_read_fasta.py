from gene_finder_YHM import read_fasta
fasta_file = "test.fasta"

description, sequence = read_fasta(fasta_file)

print(f"\nDescription: {description}")
print(f"Sequence: {sequence}")