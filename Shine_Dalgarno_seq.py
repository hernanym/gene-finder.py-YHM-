#LLM: ChatGPT(Python) Version 4o
# Prompt: Use the following code [LengthFilter_ORF.py] and 
# add a functions that looks for a ribosome binding site which is usually located 4-20bp upstream of the start codon.
# Scan upstream of the predicted start codon (e.g., -20 bp, but make it a parameter of the tool). 
# The most common ribosome binding site is the Shine-Dalgarno sequence (AGGAGG). 
# Filter all predicted ORFs based on whether they contain a Shine-Dalgarno sequence up to 20bp upstream of  the start codon.

import numpy as np
import zipfile
import os
import glob
import argparse
from Bio.Seq import Seq

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

# Mapping GCA identifiers to organism (bacteria) names
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

# Function to search for Ribosome Binding Site (Shine-Dalgarno sequence)
def has_shine_dalgarno(sequence, start_idx, upstream_max=20):
    ribosome_binding_site = "AGGAGG"
    
    # Calculate the upstream region to check (must be within sequence bounds)
    search_start = max(0, start_idx - upstream_max)
    search_end = start_idx
    
    # Get the upstream region
    upstream_region = sequence[search_start:search_end]
    
    # Check if the Shine-Dalgarno sequence is present in the upstream region
    return ribosome_binding_site in upstream_region

# Function to read FASTA file and extract sequences
def read_fasta(file_path):
    sequence = ""
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line.strip()
        return sequence
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None
    
   # Function to find ORFs using NumPy with RBS check
def find_orfs_numpy(sequence, frame=0, min_length=100, upstream_max=20): #added upstream parameter
    stop_codons = ['TAA', 'TAG', 'TGA']
    sequence = sequence[frame:]
    chars = np.array(list(sequence))
    n_codons = len(chars) // 3
    if n_codons == 0:
        return []
    codons = chars[:n_codons * 3].reshape(-1, 3)
    codon_strings = np.array([''.join(codon) for codon in codons])
    start_indices = np.where(codon_strings == 'ATG')[0]
    stop_indices = np.where(np.isin(codon_strings, stop_codons))[0]

    orfs = []
    if len(start_indices) == 0 or len(stop_indices) == 0:
        return orfs

    idx_positions = np.searchsorted(stop_indices, start_indices, side='left')

    for i, start_idx in enumerate(start_indices):
        stop_idx_pos = idx_positions[i]
        if stop_idx_pos < len(stop_indices):
            stop_idx = stop_indices[stop_idx_pos]
            if stop_idx >= start_idx:
                orf_length = stop_idx - start_idx + 1
                if orf_length >= min_length:  
                    #check for Shine-Dalgarno sequence upstream of the start codon
                    if has_shine_dalgarno(sequence, start_idx*3, upstream_max): # * 3 to map back to original sequence index
                        orf_codons = codon_strings[start_idx:stop_idx + 1]
                        orf_seq = ''.join(orf_codons)
                        orfs.append(orf_seq)
    return orfs

# Function to find ORFs in all six reading frames with min_length filter and RBS check
def find_orfs_in_six_frames(sequence, min_length=100, upstream_max=20):
    orfs = []
    frames = [0, 1, 2]
    for frame in frames:
        orfs += find_orfs_numpy(sequence, frame, min_length, upstream_max)
    reverse_complement = str(Seq(sequence).reverse_complement())
    for frame in frames:
        orfs += find_orfs_numpy(reverse_complement, frame, min_length, upstream_max)
    return orfs

# Function to translate ORFs to protein sequences
def translate_orfs(orfs):
    proteins = set()
    for orf in orfs:
        protein = Seq(orf).translate(to_stop=True)
        proteins.add(str(protein))
    return proteins


# Function to write proteins to a FASTA file
def write_proteins_to_fasta(proteins, organism_name, output_fasta_file):
    try:
        with open(output_fasta_file, 'a') as fasta_file:  # Append mode
            for i, protein in enumerate(proteins, start=1):
                header = f">{organism_name}_protein_{i}"
                fasta_file.write(f"{header}\n{protein}\n")
        print(f"Proteins successfully written to {output_fasta_file}")
    except Exception as e:
        print(f"An error occurred while writing to the file: {e}")


# Function to unzip and process genome files with RBS check and min_length filtering
def unzip_and_process_genomes(zip_path, output_fasta="SD_sequence.fasta", output_dir="extracted_genomes", min_length=100, upstream_max=20):
    try:
        # Unzip the dataset
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(output_dir)
        print(f"Successfully extracted files to {output_dir}")
    
    except Exception as e:
        print(f"Error while extracting the zip file: {e}")
        return

    # Look for all .fna files directly in the extracted directory and subdirectories
    fna_files = glob.glob(os.path.join(output_dir, '**', '*.fna'), recursive=True)

    # Filter out GCF files and keep only GCA files, ensuring each GCA ID is processed once
    processed_gca_ids = set()  # Set to track GCA IDs
    filtered_fna_files = []
    
    for fna_file in fna_files:
        file_name = os.path.basename(fna_file)
        print(f"Checking file: {file_name}")  # Debug: print the filenames being checked

        if "GCA" in file_name:
            # Use regex to capture the full GCA ID
            import re
            match = re.search(r'(GCA_\d+\.\d+)', file_name)
            if match:
                gca_id = match.group(1)  # Extract the full GCA ID
                if gca_id not in processed_gca_ids:
                    filtered_fna_files.append(fna_file)
                    processed_gca_ids.add(gca_id)
                    print(f"Added {gca_id} to the list for processing.")  # Debug: print GCA ID added

    if not filtered_fna_files:
        print("No GCA .fna files were found in the extracted dataset.")
        return

    total_files = len(filtered_fna_files)  # Get total number of files for progress reporting
    print(f"Total GCA .fna files found: {total_files}")

    # Print the full path of each file being processed
    for fna_file in filtered_fna_files:
        print(f"Found file: {fna_file}")

    for index, fna_file in enumerate(filtered_fna_files, start=1):
        file_name = os.path.basename(fna_file)
        gca_id = re.search(r'(GCA_\d+\.\d+)', file_name).group(1)  # Extract GCA ID again for safety

        # Now using the correct GCA ID to map to organism name
        organism_name = gca_to_bacteria.get(gca_id, "Unknown organism")
        print(f"\nProcessing {index} of {total_files}: {organism_name} (GCA ID: {gca_id})")

        # Step 1: Read the FASTA sequence
        sequence = read_fasta(fna_file)

        if sequence:
            # Step 2: Find ORFs in all six reading frames WITH min_length filter and Shine-Dalgarno check
            orfs = find_orfs_in_six_frames(sequence, min_length, upstream_max)  #Added upstream_max parameter

            # Step 3: Translate ORFs to protein sequences
            protein_strings = translate_orfs(orfs)

            # Step 4: Write protein sequences to the output FASTA file
            write_proteins_to_fasta(protein_strings, organism_name, output_fasta)

    print(f"Processing complete. All proteins written to {output_fasta}")

#Parse command-line arguments 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process NCBI genome dataset and filter ORFs based on Shine-Dalgarno sequence and length.")
    
    parser.add_argument("zip_path", type=str, help="Path to the NCBI dataset ZIP file")
    parser.add_argument("--min_length", type=int, default=100, help="Minimum ORF length in codons")
    parser.add_argument("--upstream_max", type=int, default=20, help="Maximum upstream distance to search for Shine-Dalgarno sequence")
    parser.add_argument("--output_fasta", type=str, default="SD_sequence.fasta", help="Output FASTA file for processed proteins")
    
    args = parser.parse_args()

    # Call the unzip_and_process_genomes function with parsed arguments
    unzip_and_process_genomes(args.zip_path, args.output_fasta, min_length=args.min_length, upstream_max=args.upstream_max)

#To Run Function on Bash as Command-Line Arguments: 
#@hernanym ➜ /workspaces/gene-finder.py-YHM- (main) $ python Shine_Dalgarno_seq.py NCBI_dataset.zip --min_length 100 --upstream_max 20 --output_fasta SD_sequence.fasta