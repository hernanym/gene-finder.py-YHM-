## Assignment Week 4 for BioE 201/230
This tool processes genome data in `.fna` format from a compressed `.zip` dataset. It scans the genomes for Open Reading Frames (ORFs), filters them by length, and checks for the presence of a Shine-Dalgarno (SD) sequence upstream of the start codon. The resulting protein sequences are then saved to a FASTA file.

## Features
- **ReadFASTA.py** - Reads a FASTA file that contains a genome using command-line arguments
- **ReadFASTA_3Frames.py** - Considers 3 Reading Frames
- **gene_finder_YHM.py** - Finds all six possible reading frames for genes and applies reverse complement.
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





