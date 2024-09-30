# LLM: ChatGPT(Python) Version 4o
# Prompt: Can you implement NumPy to process the data of my current code (NCBI_ORF_Finder.py)

import argparse
import numpy as np
import zipfile
import os
import glob
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

# Function to find ORFs using NumPy
def find_orfs_numpy(sequence, frame=0, min_length=100):
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
                    orf_codons = codon_strings[start_idx:stop_idx + 1]
                    orf_seq = ''.join(orf_codons)
                    orfs.append(orf_seq)
    return orfs

# Function to find ORFs in all six reading frames
def find_orfs_in_six_frames(sequence, min_length=100):
    orfs = []
    frames = [0, 1, 2]
    for frame in frames:
        orfs += find_orfs_numpy(sequence, frame, min_length)
    reverse_complement = str(Seq(sequence).reverse_complement())
    for frame in frames:
        orfs += find_orfs_numpy(reverse_complement, frame, min_length)
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

def unzip_and_process_genomes(zip_path, output_fasta="processed_proteins.fasta", output_dir="extracted_genomes"):
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
            # Step 2: Find ORFs in all six reading frames
            orfs = find_orfs_in_six_frames(sequence)

            # Step 3: Translate ORFs to protein sequences
            protein_strings = translate_orfs(orfs)

            # Step 4: Write protein sequences to the output FASTA file
            write_proteins_to_fasta(protein_strings, organism_name, output_fasta)

    print(f"Processing complete. All proteins written to {output_fasta}")

# Path to the NCBI dataset ZIP file and output FASTA file
#zip_file_path = "NCBI_dataset.zip"
#output_fasta_file = "processed_proteins.fasta"

# Call the function to unzip and process the genomes and write to a FASTA file
#unzip_and_process_genomes(zip_file_path, output_fasta_file)

# Main logic to parse command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process NCBI genome dataset and extract ORFs to protein sequences.")
    
    # Command-line arguments for the script
    parser.add_argument("zip_path", type=str, help="Path to the NCBI dataset ZIP file")
    parser.add_argument("--output_fasta", type=str, default="processed_proteins.fasta", help="Output FASTA file for processed proteins")
    parser.add_argument("--output_dir", type=str, default="extracted_genomes", help="Directory to extract the ZIP file")
    parser.add_argument("--min_length", type=int, default=100, help="Minimum ORF length in codons")
    
    args = parser.parse_args()
    
    # Call the unzip_and_process_genomes function with parsed arguments
    unzip_and_process_genomes(args.zip_path, args.output_fasta, args.output_dir)
    