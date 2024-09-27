#gene_finder_YHM

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


def read_fasta(file_path):

    #Reads FASTA file and returns genome sequence
    
    description = ""    # Description within FASTA fle
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

        for frame in range(3):
            possible_frames = ReadFrames(sequence, frame)
            if possible_frames: 
                print(f"\n Regions found in reading frame {frame+1}:")
                for region_data in possible_frames:
                    frame, start, end, region_sequence = region_data
                    print(f" - Region found from position {start + 1} to {end} (length: {end - start}) nucleotides")
                    print(f"    Sequence: {region_sequence[:60]}{'...' if len(region_sequence) > 60 else '' }")
            else:
                print(f"\nNo regions found in reading frame {frame + 1}.")

        return description, sequence        # Returns description (string) and complete sequence (as a single string) from function
    
    except FileNotFoundError:               # Catches error if file is not found (ie. wrong file path)
        print(f"Error: The file '{file_path}' was not found.") #prints error message
    except Exception as e:                  #Cathes unexpected errors
        print(f"An error occurred while reading the file: {e}") #prints specific error contained in e
        return None, None # If error occurs, function returns none for description and sequence
                          # meaning something went wrong.
