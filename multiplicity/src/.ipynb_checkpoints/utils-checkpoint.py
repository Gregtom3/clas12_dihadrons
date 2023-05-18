import os
import glob

def replace_elements_with_brackets(input_string, pars):
        # Replace elements in the input string with brackets
        for element in pars:
            input_string = input_string.replace(element, f'[{element}]')
        
        if "[" not in input_string:
            raise ValueError("ERROR: Unable to find variables in custom bin string")
            # If no brackets are found in the modified string, raise an error
        
        return input_string


def get_first_N_files(directory, N, glob_pattern="*"):
    # Get the list of files matching the glob pattern in the directory
    files = glob.glob(os.path.join(directory, glob_pattern))
    
    if(N>1): # Remove merged files from the pool
        files = [file for file in files if "merged" not in file]
        
    # Sort the files in alphabetical order
    sorted_files = sorted(files)

    # Get the first N files
    first_N_files = sorted_files[:N]

    return first_N_files