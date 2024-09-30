## Assignment Week 4 for BioE 201/230
This tool processes genome data in `.fna` format from a compressed `.zip` dataset. It scans the genomes for Open Reading Frames (ORFs), filters them by length, and checks for the presence of a Shine-Dalgarno (SD) sequence upstream of the start codon. The resulting protein sequences are then saved to a FASTA file.

## Features
- **ReadFASTA.py** - Takes command line parameters for a FASTA file containing a genome
- **ReadFASTA_3Frames.py** - Considers 3 Reading Frames
- **gene_finder_YHM.py** - Finds all six possible reading frames for genes and applies reverse complement.
- **OpenReadingFrames_YHM.py** - Solves Rosalind Problem (72)
- **ncbi_ORF_finder_numpy.py** - Finds all Open Reading Frames of the 14 genomes (NCBI bacteria dataset)
- **LengthFilter_ORF.py** - Takes length as a parameter and filters ORFs
- **Shine-Dalgarno Sequence Detection** - Filters ORFs based on the presence of a Shine-Dalgarno sequence upstream of the start codon.
- **FASTA Output**: Translated protein sequences are written to a FASTA file.

## Prerequisites

1. Python
2. Required Python modules:
    - `numpy`
    - `zipfile`
    - `argparse`
    - `Biopython` (`Bio.Seq`)
    - `glob`

You can install the required dependencies using `pip`:

```bash
pip install numpy biopython
```
## Code Breakdown 
### 1a. Take command line parameters for a FASTA file containing a genome 
```python
'''
LLM: ChatGPT (Python) 4o

Prompt: 
How do I write a python code that takes command-line parameters for an input file? 
The input file will consist of a FASTA file containing a genome.

'''
import argparse

def read_fasta(file_path):
    """
    Reads a FASTA file and returns the genome sequence and description.
    """
    description = ""    # Description within FASTA file
    sequence = []       # single genome sequence in file

    try:  # Starts a try-except block to handle errors that might occur while reading the file
        with open(file_path, 'r') as file:  # Opens file for reading ('r' mode)
            for line in file:  # loops through each line of the file
                line = line.strip()  # Remove leading or trailing whitespace (including newline chars)
                
                # Usually, the first line in FASTA files starts with >
                # containing the description of the sequence
                if line.startswith(">"):  # If line starts with >,
                    description = line[1:]  # Removes first character (>) and
                                            # assigns the remaining to "description" variable
                else:  # Anything else (doesn't start with >)
                    sequence.append(line)  # means it's part of the sequence and appends to sequence list
        
        sequence = ''.join(sequence)  # Concatenates all elements into a single string
        return description, sequence  # Returns description and sequence from the function
    
    except FileNotFoundError:  # Catch the error if file is not found (wrong file path)
        print(f"Error: The file '{file_path}' was not found.")
        return None, None  # If error occurs, function returns None for description and sequence
    
    except Exception as e:  # Catch unexpected errors
        print(f"An error occurred while reading the file: {e}")
        return None, None  # If error occurs, function returns None for description and sequence


# Main logic for command-line arguments
if __name__ == "__main__":
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Process a FASTA file and extract the genome sequence.")
    
    # Add an argument for the file path
    parser.add_argument("file_path", type=str, help="Path to the FASTA file")
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the read_fasta function with the file path provided by the user
    description, sequence = read_fasta(args.file_path)

    # Print the output
    if description and sequence:
        print(f"Description: {description}")
        print(f"Sequence: {sequence}")
```
### How to run using command-line parameters: 
```bash
python ReadFASTA.py test.fasta
```
### Output
```
Description: Sample Genome
Sequence: ATGCGTACGTAGCTAGTAACTGATGCGTTAACTGCGTACGATGCGTAGTAGCTGACTGAGTACGTAGTGAATGCCCGGGTACGTCGTACTGATCGTAGCTAGTAGTGCAGTACGTAGCATGACGATGCCGTAA
```

### 1b. Output any region between a start-stop codon and consider three possible reading frames. 
```python
'''
LLM: ChatGPT (Python) 4o

Prompt: 
Using ReadFASTA code, Could you make the output any region between a start ('ATG') and stop codon ('TAA", 'TAG', 'TGA')
in that FASTA file; consider three possible reading frames but ignore reverse compliments for now.

'''

import argparse

def ReadFrames (sequence, frame=0):
#Takes two arguments: 
#   - sequence = nucleotide sequence (string)
#   - frame = reading frame starting from 0 (1st nucleotide), 1 (2nd nucleotide), 2 (3rd nucleotide)

    #Finds regions within a speccific start and stop codon 
    #using 3 possible reading frames and returns a list.

    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    found_regions = [] 

    seq = sequence[frame:] # Adjusts sequence based on reading frame
    #If frame = 1, sequence will start from 2nd nucleotide, skipping 1rst nucleotide, and so forth

    #Loop over sequence in steps of 3 (codons)
    i = 0   # counter will track position in sequenc, processed in steps of 3
    while i < len(seq) - 2:  #Loop should not run if fewer than 3 nucleotides remain
        codon = seq[i:i + 3] #Extracts a codon (group of 3 nucleotides)
        if codon == start_codon: #Checks if extracted codon = to start codon

            #Found a start codon now look for a stop codon          
            for j in range(i + 3, len(seq) - 2, 3): #Starts immediately after start codon and steps through sequence 
                                                #in increments of 3 to search for stop codon, 
                                                #stops when fewer than 3 nucleotide left
                                             
                stop_codon = seq[j:j + 3]          #Extracts codon at position j to check if it is a stop codon
                if stop_codon in stop_codons:    #If it is in stop codons list
                    region_sequence = seq[i:j + 3] #Extracts the region from the start codon @ position i, to
                                                 #stop codon @ position j+3(including stop codon) and 
                                                 #stores it in variable "region_sequence"
                    found_regions.append((frame, i, j + 3, region_sequence)) #Add region to found_region list, stores 
                                                                         #reading frame, start position(i), end position(j) and region as a tuple
                    break  # Stops searching for more stop codons after first valid one is found
        i += 3 # Moves to next codon

    return found_regions # Returns list of regions found in the given reading frame

# Function to read a FASTA file
def read_fasta(file_path):
    """
    Reads a FASTA file and returns the genome sequence and description.
    """
    description = ""
    sequence = []

    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    description = line[1:]
                else:
                    sequence.append(line)
        sequence = ''.join(sequence).replace("\n", "")
        
        print(f"FASTA file description: {description}")
        print(f"Sequence length: {len(sequence)} nucleotides")
        
        return description, sequence

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return None, None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None, None


# Main function to handle reading frames and command-line arguments
def process_fasta_file(file_path):
    description, sequence = read_fasta(file_path)
    
    if description and sequence:
        for frame in range(3):
            possible_frames = ReadFrames(sequence, frame)
            if possible_frames:
                print(f"\nRegions found in reading frame {frame+1}:")
                for region_data in possible_frames:
                    frame, start, end, region_sequence = region_data
                    print(f" - Region found from position {start + 1} to {end} (length: {end - start}) nucleotides")
                    print(f"    Sequence: {region_sequence[:60]}{'...' if len(region_sequence) > 60 else ''}")
            else:
                print(f"\nNo regions found in reading frame {frame + 1}.")


if __name__ == "__main__":
    # Argument parser to handle command-line arguments
    parser = argparse.ArgumentParser(description="Process a FASTA file and find regions in different reading frames.")
    
    # Add argument for file path
    parser.add_argument("file_path", type=str, help="Path to the FASTA file")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to process the FASTA file
    process_fasta_file(args.file_path)
```
### How to run using command-line parameters: 
```bash
python ReadFASTA_3Frames.py test.fasta
```

### Output 
```
FASTA file description: Sample Genome
Sequence length: 133 nucleotides

Regions found in reading frame 1:
 - Region found from position 1 to 12 (length: 12) nucleotides
    Sequence: ATGCGTACGTAG

Regions found in reading frame 2:
 - Region found from position 22 to 30 (length: 9) nucleotides
    Sequence: ATGCGTTAA
 - Region found from position 40 to 54 (length: 15) nucleotides
    Sequence: ATGCGTAGTAGCTGA
 - Region found from position 70 to 132 (length: 63) nucleotides
    Sequence: ATGCCCGGGTACGTCGTACTGATCGTAGCTAGTAGTGCAGTACGTAGCATGACGATGCCG...
 - Region found from position 118 to 132 (length: 15) nucleotides
    Sequence: ATGACGATGCCGTAA
 - Region found from position 124 to 132 (length: 9) nucleotides
    Sequence: ATGCCGTAA

No regions found in reading frame 3.
```
### 2.  Extend tool to include the reverse complement and search all six possible reading frames

```python
LLM: ChatGPT (Python) 4o
Prompt:
I would like to extend my current code to include the reverse complement and search in all six possible reading frames for genes. I want you to use the following code for generating reverse complements of a string:  

pattern = "AAAACCCGGT"
reversePattern = pattern[::-1] 
reverseComplement = ""

for base in reversePattern:
    if base == "A":
        reverseComplement += "T"

    elif base == "T":
        reverseComplement += "A"

    elif base == "C":
        reverseComplement += "G" 

    elif base == "G":
        reverseComplement += "C"

print(reverseComplement)

How would you combine these elements?

import argparse

def ReverseComplement(sequence):
#Generates the reverse complement of a given DNA sequence

    reverse_pattern = sequence[::-1]  #Reverses the sequence
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

def ReadFrames (sequence, frame=0):
#Takes two arguments: 
#   - sequence = nucleotide sequence (string)
#   - frame = reading frame starting from 0 (1st nucleotide), 1 (2nd nucleotide), 2 (3rd nucleotide)

    #Finds regions within a speccific start and stop codon 
    #using 3 possible reading frames and returns a list.

    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    found_regions = [] 

    seq = sequence[frame:] # Adjusts sequence based on reading frame
    #If frame = 1, sequence will start from 2nd nucleotide, skipping 1rst nucleotide, and so forth

    #Loop over sequence in steps of 3 (codons)
    i = 0   # counter will track position in sequenc, processed in steps of 3
    while i < len(seq) - 2:  #Loop should not run if fewer than 3 nucleotides remain
        codon = seq[i:i + 3] #Extracts a codon (group of 3 nucleotides)
        if codon == start_codon: #Checks if extracted codon = to start codon

            #Found a start codon now look for a stop codon          
            for j in range(i + 3, len(seq) - 2, 3): #Starts immediately after start codon and steps through sequence 
                                                #in increments of 3 to search for stop codon, 
                                                #stops when fewer than 3 nucleotide left
                                             
                stop_codon = seq[j:j + 3]          #Extracts codon at position j to check if it is a stop codon
                if stop_codon in stop_codons:    #If it is in stop codons list
                    region_sequence = seq[i:j + 3] #Extracts the region from the start codon @ position i, to
                                                 #stop codon @ position j+3(including stop codon) and 
                                                 #stores it in variable "region_sequence"
                    found_regions.append((frame, i, j + 3, region_sequence)) #Add region to found_region list, stores 
                                                                         #reading frame, start position(i), end position(j) and region as a tuple
                    break  # Stops searching for more stop codons after first valid one is found
        i += 3 # Moves to next codon

    return found_regions # Returns list of regions found in the given reading frame


def read_fasta_find_regions(file_path):

    #Reads FASTA file and returns genome sequence
    
    description = ""    # Description within FASTA file
    sequence = []       # single genome sequence in file

    try:                # Starts a try-except block to handle errors
                        # that might occur while reading file

        with open(file_path, 'r') as file: # Opens file for reading ('r' mode)
            for line in file:              # loops through each line of file    
                line = line.strip()        #.strip() removes leading or triling whitespace
                                           # including newline chars
                
                #Usually first line in FASTA files start with >
                #containing the description of the sequence

                if line.startswith(">"):    # If line start with >,
                    description = line[1:]  # Removes first character (>) and 
                                            # assigns the remaining to "description" variable
                else:                       # Anything else (doesn't start with >) 
                    sequence.append(line)   # means its part of the sequence and appends to sequence list
        sequence = ''.join(sequence).replace("\n", "")       # Concatenates all elements into a single string

        print(f"FASTA file description: {description}")
        print(f"Sequence length: {len(sequence)} nucleotides")

        for frame in range(3):  #Loop iterates over the three possible reading frames on the forward strand of DNA
            found = ReadFrames(sequence, frame) #Calls ReadFrames(), "found" is a list of tuples, where each tuple represents a region found in this frame
            if found:  #Checks if any regions were found
                print(f"\n Regions found in reading frame {frame+1}:")
                for region_data in found: #inner loop goes through each region in current reading frame
                    frame, start, end, region_sequence = region_data #Unpacks tuple list (region_data) into four variables
                    print(f" - Region found from position {start + 1} to {end} (length: {end - start} nucleotides)")
                    print(f"    Sequence: {region_sequence[:60]}{'...' if len(region_sequence) > 60 else '' }") #if sequence is >60, prints the first 60 nucleotides, followed by ellipsis
            else: #if no regions were found in the current reading frame, print the following:
                print(f"\nNo regions found in reading frame {frame + 1}.")

        reverse_complement_seq = ReverseComplement(sequence) #generate reverse complement of original data sequence

        for frame in range(3): 
            #Loop similar to previous but operates on the reverse complement
            #Iterates over three reading frames: starting from 1rst nucleotide of reverse (frame=0) and so on..

            found = ReadFrames(reverse_complement_seq, frame) #Calls ReadFrames on reverse complement strand
            if found: 
                print(f"\nRegions found in reverse complement reading frame {frame + 4}:") # "frame + 4" used to represent frames as 4, 5, 6
                for region_data in found:
                    frame, start, end, region_sequence = region_data
                    print(f" - Region found from position {start + 1} to {end} (length: {end - start} nucleotides)")
                    print(f"   Sequence: {region_sequence[:60]}{'...' if len(region_sequence) > 60 else '' }")
            else: 
                print(f"\nNo regions found in reverse complement reading frame {frame + 4}.")
        
        return description, sequence
    
    except FileNotFoundError:               # Catches error if file is not found (ie. wrong file path)
        print(f"Error: The file '{file_path}' was not found.") #prints error message
    except Exception as e:                  #Cathes unexpected errors
        print(f"An error occurred while reading the file: {e}") #prints specific error contained in e
        return None, None # If error occurs, function returns none for description and sequence
                          # meaning something went wrong.

# Main function to handle command-line arguments and process the FASTA file
if __name__ == "__main__":
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Process a FASTA file and find open reading frames (ORFs).")

    # Add command-line arguments
    parser.add_argument("file_path", type=str, help="Path to the FASTA file")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function to read the FASTA file and find regions
    read_fasta_find_regions(args.file_path)
```
### How to run using command-line parameters: 
```bash
python gene_finder.py 2ndTest.fasta
```

### Output 
```
FASTA file description: TestGenomeWithReverse
Sequence length: 69 nucleotides

 Regions found in reading frame 1:
 - Region found from position 1 to 15 (length: 15 nucleotides)
    Sequence: ATGCGTAACTGCTAG

No regions found in reading frame 2.

 Regions found in reading frame 3:
 - Region found from position 34 to 54 (length: 21 nucleotides)
    Sequence: ATGCTAACGTACGGTACGTAG

No regions found in reverse complement reading frame 4.

Regions found in reverse complement reading frame 5:
 - Region found from position 1 to 30 (length: 30 nucleotides)
   Sequence: ATGCTACGTTCACTACGTACCGTACGTTAG

No regions found in reverse complement reading frame 6.
```

### 3. Rosalind - Open Reading Frame Problem
```python
LLM: ChatGPT (Python) 4o

Prompt: 
Using gene_finder_YHM.py solve the problem of identifying distinct protein strings from open reading frames (ORFs) in a given DNA sequence. Update the code so that it can translate ORFs into protein strings, consider all six reading frames and return every distinct protein string that can be derived from valid ORFs, ensuring translation stops at the first stop codon for each ORF. 

'''
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
```
### Used test_ORF.py to execute code: 
```python
from OpenReadingFrames_YHM import read_fasta_find_regions
fasta_file = "rosalind_set.fasta"

proteins = read_fasta_find_regions(fasta_file)
```
### Output
```
MDTCSAAGRVPNRPWRERQILPAAFQRSPYASSCSSAIMVPMIRR
MCRSTGTLRDFRRRSRPRGVPGVRF
MIRR
MHRTGSLGV
MSTKKRKGSFETPLAEFAAPSMVGWVPGRPRYMYPCHLLVL
MT
MGPRPGPATKLVRFGVLRHLLPSLPRYQREAVSFS
MGTMSESTF
MEVDDVGLQIVPI
MKAAKFRDVLGGGSVCDSKILRVMCAIITVTYKPKFPSESSRRQKRTPGTPRGRERLLKSLSVPVLLHIDSLADDVYGLDTVLAMPWR
MSFACIVA
M
MVYLQGIASTVSKPYTSSARESMCRSTGTLRDFRRRSRPRGVPGVRF
MEGGQTVAAVFDPPVVDLSEKRYNPLDETD
MVPMIRR
MTWIHVARPAGYPTDHGGSGKFCQRRFKGALTLLRAHLRSWFP
MERSAGVSIQST
MRLNPEQRN
MPWR
MGTMIADEHEEA
MVGWVPGRPRYMYPCHLLVL
MSISLTY
MY
MIADEHEEA
MEGAANSASGVSKEPLRFFVLICDHGSHDSPIAIGESWEPCLNQPFEQIFTFPVTSQLYVYKSNLLRIGK
MCAIITVTYKPKFPSESSRRQKRTPGTPRGRERLLKSLSVPVLLHIDSLADDVYGLDTVLAMPWR
MSISLFRIQSH
MRQTKKKKRPLVGILAVMEVDDVGLQIVPI
MYPCHLLVL
MYTVWILCLQCPGGRPWRGPLE
MSESTF
```
### 4. Find all Open Reading Frames in the 14 genomes from NCBI dataset (bacteria genomes)
```python
# LLM: ChatGPT(Python) Version 4o
# Prompt: How can I apply this code to the following genome file I downloaded from the NCBI data base, using a single command-line to find all open reading frames in all 14 genomes. I also what the code to tell me where its progressing, if there are 14 genomes in total I want it to tell me 1 of 14, 2 of 14, 3 of 14, as it progresses down the list. I want you to integrate identifiers for the files so that I know which organism is represented by which GCA file when using unzip_and_process_genomes function. Implement numpy to process the data.

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
```
### How to run using command-line parameters: 
```bash
python ncbi_ORF_finder_numpy.py NCBI_dataset.zip 
```
If you would like to specify a custom output file: 
```bash
python ncbi_ORF_finder_numpy.py NCBI_dataset.zip --output_fasta file_name.fasta
```
### Output 
```
Successfully extracted files to extracted_genomes
Checking file: GCA_000008545.1_ASM854v1_genomic.fna
Added GCA_000008545.1 to the list for processing.
Checking file: GCA_000008565.1_ASM856v1_genomic.fna
Added GCA_000008565.1 to the list for processing.
Checking file: GCA_000006825.1_ASM682v1_genomic.fna
Added GCA_000006825.1 to the list for processing.
Checking file: GCA_000027305.1_ASM2730v1_genomic.fna
Added GCA_000027305.1 to the list for processing.
Checking file: GCA_000007125.1_ASM712v1_genomic.fna
Added GCA_000007125.1 to the list for processing.
Checking file: GCA_000008745.1_ASM874v1_genomic.fna
Added GCA_000008745.1 to the list for processing.
Checking file: GCA_000008605.1_ASM860v1_genomic.fna
Added GCA_000008605.1 to the list for processing.
Checking file: GCA_000008525.1_ASM852v1_genomic.fna
Added GCA_000008525.1 to the list for processing.
Checking file: GCA_000008625.1_ASM862v1_genomic.fna
Added GCA_000008625.1 to the list for processing.
Checking file: GCA_000006745.1_ASM674v1_genomic.fna
Added GCA_000006745.1 to the list for processing.
Checking file: GCA_000008725.1_ASM872v1_genomic.fna
Added GCA_000008725.1 to the list for processing.
Checking file: GCA_000008745.1.fna
Checking file: GCA_000091085.2_ASM9108v2_genomic.fna
Added GCA_000091085.2 to the list for processing.
Checking file: GCA_000008605.1_ASM860v1_genomic.fna
Checking file: GCA_000008725.1.fna
Checking file: GCA_000008625.1_ASM862v1_genomic.fna
Checking file: GCA_000008785.1_ASM878v1_genomic.fna
Added GCA_000008785.1 to the list for processing.
Checking file: GCA_000007125.1.fna
Checking file: GCA_000006865.1_ASM686v1_genomic.fna
Added GCA_000006865.1 to the list for processing.
Checking file: GCA_000006825.1.fna
Checking file: GCA_000008745.1_ASM874v1_genomic.fna
Checking file: GCA_000008545.1.fna
Checking file: GCA_000008785.1.fna
Checking file: GCA_000008565.1.fna
Checking file: GCA_000008625.1.fna
Checking file: GCA_000008725.1_ASM872v1_genomic.fna
Checking file: GCA_000006745.1.fna
Checking file: GCA_000006825.1_ASM682v1_genomic.fna
Checking file: GCA_000007125.1_ASM712v1_genomic.fna
Checking file: GCA_000008525.1_ASM852v1_genomic.fna
Checking file: GCA_000008605.1.fna
Checking file: GCA_000008565.1_ASM856v1_genomic.fna
Checking file: GCA_000006745.1_ASM674v1_genomic.fna
Checking file: GCA_000008525.1.fna
Checking file: GCA_000091085.2.fna
Checking file: GCA_000006865.1.fna
Checking file: GCA_000027305.1_ASM2730v1_genomic.fna
Checking file: GCA_000008545.1_ASM854v1_genomic.fna
Checking file: GCA_000027305.1.fna
Checking file: GCA_000008785.1_ASM878v1_genomic.fna
Checking file: GCA_000091085.2_ASM9108v2_genomic.fna
Checking file: GCF_000027305.1_ASM2730v1_genomic.fna
Checking file: GCF_000008745.1_ASM874v1_genomic.fna
Checking file: GCA_000008545.1_ASM854v1_genomic.fna
Checking file: GCA_000008565.1_ASM856v1_genomic.fna
Checking file: GCF_000008525.1_ASM852v1_genomic.fna
Checking file: GCF_000006865.1_ASM686v1_genomic.fna
Checking file: GCA_000006825.1_ASM682v1_genomic.fna
Checking file: GCA_000027305.1_ASM2730v1_genomic.fna
Checking file: GCF_000091085.2_ASM9108v2_genomic.fna
Checking file: GCF_000008785.1_ASM878v1_genomic.fna
Checking file: GCA_000007125.1_ASM712v1_genomic.fna
Checking file: GCA_000008745.1_ASM874v1_genomic.fna
Checking file: GCA_000008605.1_ASM860v1_genomic.fna
Checking file: GCF_000006745.1_ASM674v1_genomic.fna
Checking file: GCF_000008605.1_ASM860v1_genomic.fna
Checking file: GCA_000008525.1_ASM852v1_genomic.fna
Checking file: GCF_000008625.1_ASM862v1_genomic.fna
Checking file: GCF_000008725.1_ASM872v1_genomic.fna
Checking file: GCF_000008545.1_ASM854v1_genomic.fna
Checking file: GCA_000008625.1_ASM862v1_genomic.fna
Checking file: GCA_000006745.1_ASM674v1_genomic.fna
Checking file: GCA_000008725.1_ASM872v1_genomic.fna
Checking file: GCF_000008565.1_ASM856v1_genomic.fna
Checking file: GCA_000008785.1_ASM878v1_genomic.fna
Checking file: GCF_000007125.1_ASM712v1_genomic.fna
Checking file: GCA_000091085.2_ASM9108v2_genomic.fna
Checking file: GCF_000006825.1_ASM682v1_genomic.fna
Checking file: GCA_000006865.1_ASM686v1_genomic.fna
Checking file: GCA_000006865.1_ASM686v1_genomic.fna
Total GCA .fna files found: 14
Found file: extracted_genomes/GCA_000008545.1/GCA_000008545.1_ASM854v1_genomic.fna
Found file: extracted_genomes/GCA_000008565.1/GCA_000008565.1_ASM856v1_genomic.fna
Found file: extracted_genomes/GCA_000006825.1/GCA_000006825.1_ASM682v1_genomic.fna
Found file: extracted_genomes/GCA_000027305.1/GCA_000027305.1_ASM2730v1_genomic.fna
Found file: extracted_genomes/GCA_000007125.1/GCA_000007125.1_ASM712v1_genomic.fna
Found file: extracted_genomes/GCA_000008745.1/GCA_000008745.1_ASM874v1_genomic.fna
Found file: extracted_genomes/GCA_000008605.1/GCA_000008605.1_ASM860v1_genomic.fna
Found file: extracted_genomes/GCA_000008525.1/GCA_000008525.1_ASM852v1_genomic.fna
Found file: extracted_genomes/GCA_000008625.1/GCA_000008625.1_ASM862v1_genomic.fna
Found file: extracted_genomes/GCA_000006745.1/GCA_000006745.1_ASM674v1_genomic.fna
Found file: extracted_genomes/GCA_000008725.1/GCA_000008725.1_ASM872v1_genomic.fna
Found file: extracted_genomes/NCBI_dataset/GCA_000091085.2_ASM9108v2_genomic.fna
Found file: extracted_genomes/NCBI_dataset/GCA_000008785.1_ASM878v1_genomic.fna
Found file: extracted_genomes/NCBI_dataset/GCA_000006865.1_ASM686v1_genomic.fna

Processing 1 of 14: Thermotoga maritima (GCA ID: GCA_000008545.1)
Proteins successfully written to processed_proteins.fasta

Processing 2 of 14: Deinococcus radiodurans (GCA ID: GCA_000008565.1)
Proteins successfully written to processed_proteins.fasta

Processing 3 of 14: Pasteurella multocida (GCA ID: GCA_000006825.1)
Proteins successfully written to processed_proteins.fasta

Processing 4 of 14: Haemophilus influenzae (GCA ID: GCA_000027305.1)
Proteins successfully written to processed_proteins.fasta

Processing 5 of 14: Brucella melitensis (GCA ID: GCA_000007125.1)
Proteins successfully written to processed_proteins.fasta

Processing 6 of 14: Chlamydia pneumoniae (GCA ID: GCA_000008745.1)
Proteins successfully written to processed_proteins.fasta

Processing 7 of 14: Treponema pallidum (GCA ID: GCA_000008605.1)
Proteins successfully written to processed_proteins.fasta

Processing 8 of 14: Helicobacter pylori (GCA ID: GCA_000008525.1)
Proteins successfully written to processed_proteins.fasta

Processing 9 of 14: Aquifex aeolicus (GCA ID: GCA_000008625.1)
Proteins successfully written to processed_proteins.fasta

Processing 10 of 14: Vibrio cholerae (GCA ID: GCA_000006745.1)
Proteins successfully written to processed_proteins.fasta

Processing 11 of 14: Chlamydia trachomatis (GCA ID: GCA_000008725.1)
Proteins successfully written to processed_proteins.fasta

Processing 12 of 14: Chlamydia pneumoniae (GCA ID: GCA_000091085.2)
Proteins successfully written to processed_proteins.fasta

Processing 13 of 14: Helicobacter pylori (GCA ID: GCA_000008785.1)
Proteins successfully written to processed_proteins.fasta

Processing 14 of 14: Lactococcus lactis (GCA ID: GCA_000006865.1)
Proteins successfully written to processed_proteins.fasta
Processing complete. All proteins written to processed_proteins.fasta
```
### 5. Implement a filter by length (e.g., less than 100 codons) and make the length a parameter of the tool. 
```python
LLM: # LLM: ChatGPT(Python) Version 4o
Prompt: I need to add modifications, implement a filter by length: discard short ORFs that are unlikely to be functional genes, less than 100 codons, but make the length a command-line
parameter.

# LLM: ChatGPT(Python) Version 4o
# Prompt: Use ncbi_ORF_finder_numpy.py and implement a filter by length,
# make length a parameter of the tool


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
                if orf_length >= min_length:  #USES min_length parameter for filtering
                    orf_codons = codon_strings[start_idx:stop_idx + 1]
                    orf_seq = ''.join(orf_codons)
                    orfs.append(orf_seq)
    return orfs

# Function to find ORFs in all six reading frames with min_length filter
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


#Added min_length to function parameters
def unzip_and_process_genomes(zip_path, output_fasta="filtered_proteins.fasta", output_dir="extracted_genomes", min_length=100):
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
            # Step 2: Find ORFs in all six reading frames WITH min_length filter
            orfs = find_orfs_in_six_frames(sequence, min_length)  

            # Step 3: Translate ORFs to protein sequences
            protein_strings = translate_orfs(orfs)

            # Step 4: Write protein sequences to the output FASTA file
            write_proteins_to_fasta(protein_strings, organism_name, output_fasta)

    print(f"Processing complete. All proteins written to {output_fasta}")

# Path to the NCBI dataset ZIP file and output FASTA file
zip_file_path = "NCBI_dataset.zip"
output_fasta_file = "filtered_proteins.fasta"

# Call the function to unzip and process the genomes and write to a FASTA file
unzip_and_process_genomes(zip_file_path, output_fasta_file)
```
### To run the script with a different minimum ORF length, 
Simply adjust the min_length parameter.
### To Run on BASH: 
```python 
python LengthFilter_ORF.py --min_length 100
```
### Output is a FASTA file
```
... same
.. as 
. code above

Proteins successfully written to filtered_proteins.fasta
Processing complete. All proteins written to filtered_proteins.fasta
```
### 6. Implement a filter by length (e.g., less than 100 codons) and make the length a parameter of the tool. 
```python

# LLM: ChatGPT(Python) Version 4o
# Prompt: Use the following code {LLengthFilter_ORF.py} add a function that looks for a ribosome binding site which is usually located 4-20bp upstream of the start codon. Scan upstream of the predicted start codon (e.g., -20 bp) and make this a parameter. Filter all predicted ORFs based on whether they contain a Shine-Dalgarno sequence up to 20bp upstream of the start codon.

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
```
### To Run Function on Bash as Command-Line Argument:
```bash
python Shine_Dalgarno_seq.py NCBI_dataset.zip --min_length 100 --upstream_max 20 --output_fasta SD_sequence.fasta
```
### Output is a FASTA file 
```
... same
.. as 
. code above

Proteins successfully written to SD_sequence.fasta
Processing complete. All proteins written to SD_sequence.fasta
```


