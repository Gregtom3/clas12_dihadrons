import os
import glob
import re
def retitle(var):
    
    VAR = var.replace("rec_","").replace("gen_","")
    
    if VAR == "z":
        return "z"
    elif VAR == "Mh":
        return "M_{h} [GeV]"
    elif VAR == "x":
        return "x"
    elif VAR == "Q2":
        return "Q^{2} [GeV^{2}]"
    elif VAR == "y":
        return "y"
    elif VAR == "pTtot":
        return "p_{T} [GeV]"
    elif VAR == "xF":
        return "x-Feynman"
    elif VAR == "xF1":
        return "x-Feynman(h_{1})"
    elif VAR == "xF2":
        return "x-Feynman(h_{2})"
    elif VAR == "z1":
        return "z(h_{1})"
    elif VAR == "z2":
        return "z(h_{2})"
    elif VAR == "phi_h":
        return "\phi_{h}"
    elif VAR == "phi_R0":
        return "\phi_{R}"
    elif VAR == "th":
        return "\theta_{COM}"
    else:
        return var

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

def reformat_rect_bins(bin_edges,bin_names):
    # Ensure the bin edges and names are of same length
    if len(bin_edges) != len(bin_names):
        raise ValueError("Number of bin edge lists and bin names must be equal")

    result_strings = []

    # Create a list of ranges for each bin
    bin_ranges = [list(zip(edges[:-1], edges[1:])) for edges in bin_edges]

    # Generate a Cartesian product of the ranges to get all combinations
    from itertools import product
    for combination in product(*bin_ranges):
        bin_string = ''
        for name, (lower, upper) in zip(bin_names, combination):
            bin_string += f'({name}>{lower})&({name}<{upper})&'
        # Remove trailing '&' and add to the result list
        result_strings.append(bin_string[:-1])

    return result_strings

def split_into_fours(lst):
    max_groups = max(len(s.split('&')) for s in lst)
    ret = [[] for _ in range((max_groups + 3) // 4)]  # create a list of lists based on max_groups
    for string in lst:
        groups = string.split('&')
        for i in range(0, len(groups), 4):
            ret[i // 4].append("&".join(groups[i:i+4]))
    return ret

def remove_duplicates(lst):
    unique_list = []
    for item in lst:
        if item not in unique_list:
            unique_list.append(item)
    return unique_list

def split_into_twos(lst):
    return [lst[i:i+2] for i in range(0, len(lst), 2)]

def extract_number_with_decimal(string):
    pattern = r"\d+(\.\d+)?"
    match = re.search(pattern, string)
    if match:
        return float(match.group())
    else:
        return None
    
    
def calculate_iou(hist1, hist2):
    # Get the bin contents for both histograms
    bins1 = [hist1.GetBinContent(i) for i in range(1, hist1.GetNbinsX() + 1)]
    bins2 = [hist2.GetBinContent(i) for i in range(1, hist2.GetNbinsX() + 1)]

    # Calculate the intersection and union
    intersection = sum(min(bins1[i], bins2[i]) for i in range(len(bins1)))
    union = sum(max(bins1[i], bins2[i]) for i in range(len(bins1)))

    # Calculate the IoU
    try:
        iou = intersection / union
    except:
        iou = 0
    return iou