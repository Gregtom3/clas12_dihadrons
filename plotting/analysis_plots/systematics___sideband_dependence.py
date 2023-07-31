######################################################################
import uproot
import os
import argparse
import matplotlib.pyplot as plt
from tools___asym import *
from tools___io import *
from tools___etc import *
import copy

# Create parser
parser = argparse.ArgumentParser(
                prog='systematics___sideband_dependence.py',
                description='Plot the asymmetries for dihadrons in different Mgg regions')
# Add the optional recreate flag
parser.add_argument('--twoh', action='store_true', help='Use 9 modulations (2h)')
args = parser.parse_args()

# Enable the LaTeX rendering backend
plt.rcParams["text.usetex"] = True

#####################################################################
# Program: systematics___sideband_dependence.py
# Author:  Gregory Matousek
# Date:    7/17/2023
#
# Purpose: Calculate+Plot asymmetries in different regions of M_gg
# Requirements: Run `scripts/calc_dihadron_asym_precut_many_sideband.sh`
#####################################################################

def main():
    def make_plot(par,ax,outtext=True):
        # Create figure
        #fig,ax = plt.subplots(1,1,dpi=150,figsize=(4,4))

        # Loop over each sideband experiment
        for i,sideband_dir in enumerate(sideband_dirs):
            # Get the Mgg_min and Mgg_max from the dir
            Mgg_string = sideband_dir.split("/")[-1]
            Mgg_min = float(Mgg_string.split("_")[0])
            Mgg_max = float(Mgg_string.split("_")[1])
            Mgg_mid = 0.5*(Mgg_min+Mgg_max)
            Mgg_width = 0.5*(Mgg_max-Mgg_min)

            # Open the Result root file using Uproot
            try:
                result = uproot.open(sideband_dir + f"/outObsBins_sdbnd/{file_name}")
                resultTree = result[ttree_name]

                # Get the "A" and "A_err" branches using Uproot
                branch = resultTree[par]
                err_branch = resultTree[f"{par}_err"]

                # Convert the branches to numpy arrays
                value = branch.array()[0]
                err = err_branch.array()[0]

                ax.errorbar(Mgg_mid,value,xerr=Mgg_width,yerr=err,fmt="ko",capsize=3)
            except:
                continue

        # Create a twin axis for histogram
        axHist = ax.twinx()

        # Make a faint histogram of Mgamma
        pion_pair = pull_pair_from_string(file_path)
        if pion_pair=="piplus_pi0":
            u = uproot.open("/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/pipi0_paper_RGA_only/volatile/data/piplus_pi0/Fall2018Spring2019_RGA_inbending_merged_cuts.root")
        else:
            u = uproot.open("/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/pipi0_paper_RGA_only/volatile/data/piminus_pi0/Fall2018_RGA_outbending_merged_cuts.root")
        u = u["dihadron_cuts"]
        Mgg_data = u["M2"].array(library="np")
        #noML_data =  u["isGoodEventWithoutML"].array(library="np")
        #Mgg_data = Mgg_data[noML_data==1]
        # Determine maximum height of histogram
        n, bins, patches = axHist.hist(Mgg_data, bins=100,range=(0.01,0.8), alpha=0.5, density=True)
        hist_height = n.max()

        # Determine the ratio of the desired histogram height (in relation to the main y-axis) and the actual height
        scale_factor = ax.get_ylim()[1] / hist_height
        axHist.axis('off')  # Hide the secondary axis
        # Set symmetric y-axis limits
        #y_max = max(abs(np.array(ax.get_ylim())))
        y_max=0.045
        ax.set_ylim(-y_max, y_max)
        # Add y-axis grid
        ax.yaxis.grid(True)
        # Set x-axis label
        ax.set_xlabel("$M_{\gamma\gamma}[GeV]$",fontsize=15)
        # Set y-axis label
        ax.set_ylabel("$A_{LU}^{"+mod_data[par]+"}$",fontsize=15)
        # Add darker line at y=0
        ax.axhline(0, color='black', linewidth=1)
        # Add text
        version = pull_version_from_string(file_path)
        version = printable_version(version)
        if outtext:
            ax.text(0.7, 1.03, f"${latex_pion_pair(pion_pair)}$ {version}", fontsize=12, color='black', ha='center', transform=plt.gca().transAxes)





    # Loop over all sideband subdirectories
    for file_path,file_name,ttree_name in zip(io.file_paths,
                                              io.file_names,
                                              io.ttree_names):
        
        # Pull list of sideband runs
        sideband_dirs = sorted([file_path+d for d in os.listdir(file_path)])
        
        # Save modulation data
        if args.twoh:
            mod_data = get_2h_modulations(2)
        else:
            mod_data = get_modulations(2)

        # List of all pars
        pars = mod_data.keys()
        # Define number of rows and columns for the subplot grid
        n_cols = 3
        n_rows = int(np.ceil(len(pars)/n_cols))  # assuming we have enough pars


        # Create a figure for the combined plots
        bigfig, axs = plt.subplots(n_rows, n_cols, dpi=150, figsize=(4*n_cols, 4*n_rows))

        # Ensure axs is always a 2D array
        if n_rows == 1 and n_cols == 1:
            axs = np.array([[axs]])
        elif n_rows == 1 or n_cols == 1:
            axs = axs.reshape(n_rows, n_cols)
        
        for idx,par in enumerate(pars):
            fig, ax = plt.subplots(1,1,dpi=150,figsize=(4,4))
            make_plot(par,ax)
            make_plot(par,axs[idx//n_cols, idx%n_cols],outtext=False)
            pion_pair = pull_pair_from_string(file_path)
            fig.savefig(io.out+"/"+pion_pair+f"/{par}.png",bbox_inches='tight')
            print("Figure saved:",io.out+"/"+pion_pair+f"/{par}.png")
            
        # Manually adjust layout
        bigfig.subplots_adjust(hspace = 0.5, wspace = 0.5)
        for ax in axs.flatten():
            if not ax.lines:
                bigfig.delaxes(ax)
                
        # Add text
        version = pull_version_from_string(file_path)
        version = printable_version(version)
        #bigfig.text(0.75, 0.15, f"${latex_pion_pair(pion_pair)}$ {version}", fontsize=20, color='black', ha='center')#, transform=plt.gca().transAxes)  
        bigfig.savefig(io.out+"/"+pion_pair+f"/aggregate.png",bbox_inches='tight')
        print("Aggregate figure saved:", io.out+"/"+pion_pair+f"/aggregate.png")
        
if __name__ == "__main__":
    
    global io
    
    io = DataIO(PROJECT_NAME="pipi0_paper_RGA_only",
                file_paths=["/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/pipi0_paper_RGA_only/systematics/many_sideband_2h_ML/Fall2018Spring2019_RGA_inbending/precut/piplus_pi0/", "/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/pipi0_paper_RGA_only/systematics/many_sideband_2h_ML/Fall2018_RGA_outbending/precut/piminus_pi0/"],
                file_names=["ResultsHSMinuit2.root"]*2,
                ttree_names=["ResultTree"]*2,
                output_dir="sideband_dependence_2h_ML")
    
    io.create_multiple_pion_pair_dirs()
    


    main()
