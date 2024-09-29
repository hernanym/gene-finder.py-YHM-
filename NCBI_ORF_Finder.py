import zipfile #module provides functions and classes for working with ZIP archives
import os  #module in Python that provides a way of interacting with operating system
import glob #module is used to find files or directories that match a specified pattern 

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

#Mapping GCA identifiers to organism (bacteria) names 
gca_to_bacteria = {
    "GCA_000006745.1": "Vibrio cholerae",
    "GCA_000006825.1": "Pasteurella multocida",
    "GCA_000006865.1": "Lactococcus lactis",
    "GCA_000007125.1": "Brucella melitensis",
    "GCA_000008525.1": "Helicobacter pylori",
    "GCA_000008545.1": "Thermotoga maritima",
    "GCA_000008565.1": "Deinococcus radiodurans",
    "GCA_000008605.1": "Treponema pallidum",
    "GCA_000008625.1": "Aquifex aeolicus",
    "GCA_000008725.1": "Chlamydia trachomatis",
    "GCA_000008745.1": "Chlamydia pneumoniae",
    "GCA_000008785.1": "Helicobacter pylori",
    "GCA_000091085.2": "Chlamydia pneumoniae",
    "GCA_000027305.1": "Haemophilus influenzae"
}
      


def read_fasta_find_regions(file_path, organism_name):
    print(f"Reading and processing {file_path} for organism: {organism_name}...")

    description = ""  # Description from the FASTA file header
    sequence = []  # List to store the DNA sequence

    try:
        # Open the FASTA file for reading
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace (including newlines)
                if line.startswith(">"):  # If the line starts with '>', it's the description
                    description = line[1:]  # Remove the '>' and store the description
                else:
                    sequence.append(line)  # Append DNA sequence lines to the list

        # Combine all sequence lines into one string
        sequence = ''.join(sequence).replace("\n", "")

        # Print the description here (or store it with other data)
        print(f"FASTA Description: {description}")

        proteins = set()  # Initialize an empty set to store distinct proteins

        # Search ORFs in the forward strand (reading frames 1, 2, and 3)
        for frame in range(3):
            proteins.update(ReadFrames(sequence, frame))

        # Generate reverse complement of the sequence to search reverse reading frames
        reverse_complement_seq = ReverseComplement(sequence)

        # Search ORFs in the reverse complement strand (reading frames 4, 5, and 6)
        for frame in range(3):
            proteins.update(ReadFrames(reverse_complement_seq, frame))

        # Print each distinct protein string found in the sequence along with organism name
        for protein in proteins:
            print(f"Organism: {organism_name}, Protein: {protein}")

        return proteins  # Return the set of distinct proteins
    
    except FileNotFoundError:
        # Handle the case where the file is not found
        print(f"Error: The file '{file_path}' was not found.")
        return None
    except Exception as e:
        # Handle any other exceptions that may occur during file reading
        print(f"An error occurred while reading the file: {e}")
        return None


def unzip_and_process_genomes(zip_path, output_dir="extracted_genomes"):
    #Unzips the dataset and applies ORF finding to all genome files within the unzipped directory.
      
    # Step 1: Unzip the dataset
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Extract all the files in the ZIP archive into the output directory
            zip_ref.extractall(output_dir)
        print(f"Successfully extracted files to {output_dir}")
    
    except Exception as e:
        # Handle any errors during the unzipping process
        print(f"Error while extracting the zip file: {e}")
        return

    # Step 2: Find all genome files in the extracted directory (assuming they have .fna extension)
    fna_files = glob.glob(os.path.join(output_dir, '**', '*.fna'), recursive=True)
    
    if not fna_files:
        # If no .fna files are found, print a message and exit
        print("No FNA files were found in the extracted dataset.")
        return

    total_files = len(fna_files)  # Get total number of files for progress reporting
    all_proteins = []   #List to store all proteins from all genomes with with their associated info

    # Step 3: Loop through each genome file and apply the ORF finder
    for index, fna_file in enumerate(fna_files, start=1):
        # Extract the GCA identifier from the file name
        gca_id = os.path.basename(fna_file).split("_")[0]

        # Look up the organism name using the GCA identifier
        organism_name = gca_to_bacteria.get(gca_id, "Unknown organism")

        # Print progress with the organism name
        print(f"\nProcessing {index} of {total_files}: {organism_name} (GCA ID: {gca_id})")

        proteins_with_source = read_fasta_find_regions(fna_file, organism_name)

        # Add the results to the overall list of proteins
        all_proteins.extend(proteins_with_source)

    if all_proteins: 
        print("\nProtein sequences found:")
        for protein, file_path, organism_name in all_proteins:
            print(f"Protein: {protein}\nFile: {file_path}\nOrganism: {organism_name}\n")
    else: 
        print("No protein sequences were found in any genome files.")

        # Apply the ORF finding function to the current genome file
        #read_fasta_find_regions(fna_file)

# Path to the NCBI dataset ZIP file (this should be updated to the correct path in Codespaces)
zip_file_path = "ncbidataset.zip"  # This is where your uploaded ZIP file is stored

# Call the function to unzip and process the genomes
unzip_and_process_genomes(zip_file_path)
