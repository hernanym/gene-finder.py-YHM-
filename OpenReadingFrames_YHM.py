# Codon table for DNA to protein translation
codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'STOP', 'TAG': 'STOP', 'TGC': 'C', 'TGT': 'C', 'TGA': 'STOP', 'TGG': 'W',
}

def ReverseComplement(sequence):
    #Generates the reverse complement of a given DNA sequence

    reverse_pattern = sequence[::-1]  # Reverses the sequence
    reverse_complement = "" 

    for base in reverse_pattern:
        if base == "A":
            reverse_complement += "T"
        elif base == "T":
            reverse_complement += "A"
        elif base == "C":
            reverse_complement += "G"
        elif base == "G":
            reverse_complement += "C"
    
    return reverse_complement

def translate_orf_to_protein(orf_sequence):
    # Translates an ORF (DNA sequence) into a protein sequence.

    protein = [] #List to store the translated amino acids

    #Iterates through DNA sequence in chunks of codons
    for i in range(0, len(orf_sequence) - 2, 3):
        codon = orf_sequence[i:i + 3]  #Extract a codon
        amino_acid = codon_table.get(codon, "") #Get corresponding amino acid from codon table

        if amino_acid == "STOP": #If a stop codon is encountered, stop translation
            break

        if amino_acid:          #If it's a valid codon, append the amino acid to the protein list
            protein.append(amino_acid)

    return "".join(protein)     #Return the protein sequence as a string

def ReadFrames(sequence, frame=0):
    #Find ORFs in a specific reading frame and return protein strings

    start_codon = 'ATG'   #start codon for protein translation
    stop_codons = {'TAA', 'TAG', 'TGA'} #set of possible stop codons
    found_proteins = set()  # Set to store distinct protein strings

    seq = sequence[frame:]  # Adjust sequence based on reading frame
    i = 0

    while i < len(seq) - 2: #Loop through sequence in codon steps
        codon = seq[i:i + 3]

        if codon == start_codon: #Extracts codon
            # Found a start codon, now look for a stop codon
            for j in range(i + 3, len(seq) - 2, 3): #Start searching after start codon
                stop_codon = seq[j:j + 3]  #Extract the next codon
                if stop_codon in stop_codons: #If codon is stop codon
                    orf_sequence = seq[i:j + 3] #Extract the entire ORF (from start-stop codon)
                    protein = translate_orf_to_protein(orf_sequence)
                    if protein:
                        found_proteins.add(protein)  # Add the protein to the set 
                    break #Stop searching after finding the first valid stop codon

        i += 3  # Move to the next codon

    return found_proteins #Return the set of proteins found in this reading frame

def read_fasta_find_regions(file_path):
    #Reads a FASTA file and returns distinct protein strings from all ORFs

    description = ""    # Description within FASTA file
    sequence = []       # DNA sequence from file

    try:  # Try-except block to handle errors during file reading
        with open(file_path, 'r') as file:  # Open the file for reading
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace, including newlines
                if line.startswith(">"):
                    description = line[1:]  # Remove the '>' and store the description
                else:
                    sequence.append(line)  # Append sequence lines
        sequence = ''.join(sequence).replace("\n", "")  # Combine the sequence into a single string

        proteins = set()  # Store distinct protein strings

        # Search ORFs in the forward strand (reading frames 1, 2, 3)
        for frame in range(3):
            proteins.update(ReadFrames(sequence, frame))

        # Generate reverse complement of the sequence
        reverse_complement_seq = ReverseComplement(sequence)

        # Search ORFs in the reverse complement strand (reading frames 4, 5, 6)
        for frame in range(3):
            proteins.update(ReadFrames(reverse_complement_seq, frame))

        # Print each distinct protein string
        for protein in proteins:
            print(protein)

        return proteins  # Return the set of distinct proteins
    
    except FileNotFoundError:
        #Handles errors when file not  found
        print(f"Error: The file '{file_path}' was not found.")  # Print file not found error
    except Exception as e:
        #Handles any other error that may occur while reading file
        print(f"An error occurred while reading the file: {e}")  # Print any other error encountered
        return None, None
