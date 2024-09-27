#gene_finder_YHM

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
        sequence = ''.join(sequence)        # Concatenates all elements into a single string
        return description, sequence        # Returns description (string) and complete sequence (as a single string) from function
    
    except FileNotFoundError:               # Catches error if file is not found (ie. wrong file path)
        print(f"Error: The file '{file_path}' was not found.") #prints error message
    except Exception as e:                  #Cathes unexpected errors
        print(f"An error occurred while reading the file: {e}") #prints specific error contained in e
        return None, None # If error occurs, function returns none for description and sequence
                          # meaning something went wrong.
