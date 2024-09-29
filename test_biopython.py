from Bio.Seq import Seq

# Create a DNA sequence
dna_seq = Seq("AGTACACTGGT")

# Get the reverse complement of the sequence
reverse_complement = dna_seq.reverse_complement()

print(reverse_complement)
