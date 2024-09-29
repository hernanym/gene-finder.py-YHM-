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