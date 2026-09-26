import argparse
import yaml
import ROOT
import os
from array import array
import numpy as np
import scipy
import copy

# Function to generate logarithmically spaced bins
def logspace_bins(start, stop, num):
    return np.logspace(np.log10(start), np.log10(stop), num)

def get_4d_hists(infile,ttree,xQ2_boundaries, zpT_boundaries, var1, var2, var1_bins, var2_bins, var3, var3_bins, var4=None, var4_bins=None):
    
    # create RDataFrame and select events that satisfy the cut conditions
    df = ROOT.RDataFrame(ttree, infile)
    df_cut=df
    
    # create main histogram
    main_hist = df_cut.Histo2D(("h2d_main_{0}_{1}".format(var1, var2), "{0}-{1} distribution".format(var1, var2), len(var1_bins) - 1, array('d', var1_bins), len(var2_bins) - 1, array('d', var2_bins)), var1, var2)

    
    var1_var2_hists = []
    var3_var4_hists = []
    var3_var4_subhists = []
    
    tmp_df_cuts=[]
    for idx, b in enumerate(xQ2_boundaries):
        print(idx + 1, "of", len(xQ2_boundaries))
        tmp_df_cut=df_cut.Filter("&&".join(b))
        tmp_df_cuts.append(tmp_df_cut)
        
        # create histogram for var1 and var2
        h2d_var1_var2 = tmp_df_cuts[-1].Histo2D(("h2d_{0}_{1}_{2}".format(var1, var2, idx), "{0}-{1} distribution".format(var1, var2), len(var1_bins) - 1, array('d', var1_bins), len(var2_bins) - 1, array('d', var2_bins)), var1, var2)

        var1_var2_hists.append(h2d_var1_var2)

        h2d_var3_var4 = tmp_df_cuts[-1].Histo2D(("h2d_{0}_{1}_{2}".format(var3, var4, idx), "{0}-{1} distribution".format(var3, var4), len(var3_bins) - 1, array('d', var3_bins), len(var4_bins) - 1, array('d', var4_bins)), var3, var4)

            
        var3_var4_hists.append(h2d_var3_var4)
        
        var3_var4_subhists_arr=[]
        
        for idxidx,bb in enumerate(zpT_boundaries):
            tmptmp_df_cut = tmp_df_cut.Filter("&&".join(bb))
            tmp_df_cuts.append(tmptmp_df_cut)
            h2d_var3_var4_sub = tmp_df_cuts[-1].Histo2D(("h2d_{0}_{1}_{2}_sub".format(var3, var4, idx), "{0}-{1} distribution".format(var3, var4), len(var3_bins) - 1, array('d', var3_bins), len(var4_bins) - 1, array('d', var4_bins)), var3, var4)

            var3_var4_subhists_arr.append(h2d_var3_var4_sub)
    
        var3_var4_subhists.append(var3_var4_subhists_arr)

    
    return main_hist.GetValue(), [h.GetValue() for h in var1_var2_hists], [h.GetValue() for h in var3_var4_hists], [[h.GetValue() for h in hlist] for hlist in var3_var4_subhists]



def get_2d_hists(infile,ttree,xQ2_boundaries, var1, var2, var1_bins, var2_bins):
    
    # create RDataFrame and select events that satisfy the cut conditions
    df = ROOT.RDataFrame(ttree, infile)
    df_cut=df
    
    # create main histogram
    main_hist = df_cut.Histo2D(("h2d_main_{0}_{1}".format(var1, var2), "{0}-{1} distribution".format(var1, var2), len(var1_bins) - 1, array('d', var1_bins), len(var2_bins) - 1, array('d', var2_bins)), var1, var2)

    
    var1_var2_hists = []
    
    tmp_df_cuts=[]
    for idx, b in enumerate(xQ2_boundaries):
        print(idx + 1, "of", len(xQ2_boundaries))
        tmp_df_cut=df_cut.Filter("&&".join(b))
        tmp_df_cuts.append(tmp_df_cut)
        
        # create histogram for var1 and var2
        h2d_var1_var2 = tmp_df_cuts[-1].Histo2D(("h2d_{0}_{1}_{2}".format(var1, var2, idx), "{0}-{1} distribution".format(var1, var2), len(var1_bins) - 1, array('d', var1_bins), len(var2_bins) - 1, array('d', var2_bins)), var1, var2)

        var1_var2_hists.append(h2d_var1_var2)
    
    return main_hist.GetValue(), [h.GetValue() for h in var1_var2_hists]

def generate_x_Q2_boundaries(x_bins, add_middle_cut):
    input_boundaries = []
    
    for i in range(len(x_bins)):
        x_min=x_bins[i][0]
        x_max=x_bins[i][1]
        input_boundary = ["x<={0}".format(x_max), "x>={0}".format(x_min)]
        if add_middle_cut[i] == -1:
            input_boundary.append("Q2<{0}".format(middle_Q2_cut))
        elif add_middle_cut[i] == 1:
            input_boundary.append("Q2>{0}".format(middle_Q2_cut))

        input_boundary.append(high_Q2_cut)
        input_boundary.append(low_Q2_cut)
        input_boundaries.append(input_boundary)
        
    return input_boundaries

def generate_z_pt_boundaries(z_edges, pT_edges):
    boundaries = []
    nzbins = len(z_edges) - 1
    npTbins = len(pT_edges) - 1

    for iz in range(nzbins):
        for ip in range(npTbins):
            boundaries.append(["z>={0}".format(z_edges[iz]), "z<{0}".format(z_edges[iz+1]), "pTtot>={0}".format(pT_edges[ip]), "pTtot<{0}".format(pT_edges[ip+1])])


    return boundaries

def create_nested_4d_dict(x_arr, z_arr, pT_arr, xQ2_hists, zpT_subhists):
    # Initialize an empty dictionary to store the nested counts
    nested_dict = {}

    # Loop over each value in x_arr
    for i in range(len(x_arr)):
        # Initialize an empty dictionary for this value of i
        name_i = "x_Q2_bin_{}".format(i)
        x = np.round(xQ2_hists[i].GetMean(1),4)
        Q2 = np.round(xQ2_hists[i].GetMean(2),4)
        
        nested_dict[name_i] = {"x":float(x),"Q2":float(Q2)}

        # Loop over each value in z_arr, except for the last value
        for j in range(len(z_arr)-1):
            # Initialize an empty dictionary for this value of j
            name_j = "z_bin_{}".format(j)
            z = z_arr[j]
            nested_dict[name_i][name_j] = {"zmin":float(z_arr[j]),"zmax":float(z_arr[j+1]),"z":float(0.5*(z_arr[j+1]+z_arr[j]))}

            # Loop over each value in pT_arr, except for the last value
            for k in range(len(pT_arr)-1):
                # Initialize an empty dictionary for this value of k
                name_k = "pT_bin_{}".format(k)
                pT = np.round((pT_arr[k] + pT_arr[k+1])/2,3)

                # Calculate the bin number for this (z,pT) pair
                zpT_binnum = j*(len(pT_arr)-1)+k

                # Get the number of counts for this (z,pT) pair
                counts = zpT_subhists[i][zpT_binnum].GetEntries()

                # Store the counts and pT value in the nested dictionary
                nested_dict[name_i][name_j][name_k] = {"Counts":counts,"Error":float(np.sqrt(counts)),"pT":float(pT),"pTmin":float(pT_arr[k]),"pTmax":float(pT_arr[k+1])}

    # Return the completed nested dictionary
    return nested_dict



def create_nested_2d_dict(x_arr, xQ2_hists):
    # Initialize an empty dictionary to store the nested counts
    nested_dict = {}

    # Loop over each value in x_arr
    for i in range(len(x_arr)):
        # Initialize an empty dictionary for this value of i
        name_i = "x_Q2_bin_{}".format(i)
        x = np.round(xQ2_hists[i].GetMean(1),4)
        Q2 = np.round(xQ2_hists[i].GetMean(2),4)
        
        # Get the number of counts for this (x,Q2) bin
        counts = xQ2_hists[i].GetEntries()

        # Store the counts
        nested_dict[name_i] = {"x":float(x),"Q2":float(Q2), "Counts":counts,"Error":float(np.sqrt(counts))}

    # Return the completed nested dictionary
    return nested_dict



# Define the global boundary expressions
low_Q2_cut = "Q2>1.4144 + -5.4708 * x + 40.5357 * x*x + -40.0208 * x*x*x + 29.2121 * x*x*x*x"
middle_Q2_cut = "0.6361935324019532 + 5.961630973508846*x + 12.028695097029118*x*x"
high_Q2_cut = "Q2<x*17"

# Main function
def main():
    
    # Define the argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Input .root file")
    parser.add_argument("output_file", help="Output .yaml file")
    parser.add_argument("version", help="[dihadron] or [dis]")
    
    # Parse the command-line arguments
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    version = args.version
    # Check that the input file exists
    if not os.path.isfile(args.input_file):
        print("Error: input file {} does not exist".format(args.input_file))
        return
    
    # Check that the output file ends in .yaml
    if not args.output_file.endswith(".yaml"):
        print("Error: output file {} must end in .yaml".format(args.output_file))
        return
    
    # Check that the version is valid
    if args.version not in ["dihadron", "dis"]:
        print("Error: version must be dihadron or dis")
        return
    

    ttree  = ("dihadron_cuts" if version=="dihadron" else "DIS")
    
    # Define the boundary expressions
    #x_arr = [[0.07,0.12],[0.12,0.2],[0.2,0.275],[0.275,0.42],[0.42,1], [0.12,0.15],[0.15,0.22],[0.22,0.29],[0.29,0.42]]
    #Q2_arr = [0,-1,-1,-1,0,1,1,1,1] # Same length as x arr. 0 --> dont use center line, -1 below center line, +1 above center line
    x_arr = [[0.07,0.12], [0.12,0.15], [0.12,0.2],[0.15,0.22],[0.2,0.275],[0.22,0.29],[0.275,0.42],[0.29,0.42],[0.42,1]]
    Q2_arr = [0,1,-1,1,-1,1,-1,1,0] # Same length as x arr. 0 --> dont use center line, -1 below center line, +1 above center line
    
    z_arr = [0.3,0.38,0.5,0.6,0.7,0.9]
    pT_arr = [0,0.2,0.4,0.6,0.8,1,1.5]

    xQ2_boundaries = generate_x_Q2_boundaries(x_arr, Q2_arr)
    zpT_boundaries = generate_z_pt_boundaries(z_arr,pT_arr)

    x_bins = logspace_bins(5e-2, 1, 100)
    Q2_bins = logspace_bins(1, 20, 100)
    z_bins = np.linspace(0, 1, 100)
    pt_bins = np.linspace(-0.3, 2, 100)

    print("Generating Histograms from TFile")
    if(ttree=="dihadron_cuts"): # 4d
        xQ2_hist_main, xQ2_hists, zpT_hists, zpT_subhists = get_4d_hists(input_file, ttree, xQ2_boundaries, zpT_boundaries,"x", "Q2", x_bins, Q2_bins, "z", z_bins, "pTtot", pt_bins)
        print("Making YAML dictionary")
        my_nested_dict = create_nested_4d_dict(x_arr, z_arr, pT_arr, xQ2_hists,  zpT_subhists)
    else: # 2d
        xQ2_hist_main, xQ2_hists = get_2d_hists(input_file, ttree, xQ2_boundaries,"x", "Q2", x_bins, Q2_bins)
        print("Making YAML dictionary")
        my_nested_dict = create_nested_2d_dict(x_arr, xQ2_hists)
    
    # Get the parent directory path
    parent_dir = os.path.dirname(output_file)

    # Make the parent directory if it doesn't exist
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)
        
    # Write the dictionary to the output file in YAML format
    with open(output_file, "w") as f:
        yaml.dump(my_nested_dict, f)
    
    print("Done")
    
if __name__ == '__main__':
    main()