#######################################################################
import itertools
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
})
plt.style.use("science")
import matplotlib as mpl
import argparse
from tools___asym import *



# Create parser
parser = argparse.ArgumentParser(
                prog='asymmetry__seven_modulations.py',
                description='Create asymmetry plots for the dihadron analysis (integrate over dihadron theta, so not partial waves)')


#####################################################################
# Program: asymmetry___seven_modulations.py
# Author:  Gregory Matousek
# Date:    7/17/2023
#
# Purpose: Plot asymmetries of 7 phi_h/phi_R sinusoidal modulations
# Requirements: Run `scripts/calc_dihadron_asym_precut.sh`
#####################################################################

def main():
    
    # Create new YAML file
    if args.create == True:
        print("Not running `create_asym_yaml` just in case")
        #create_asym_yaml(io.PROJECT_DIR,io.PROJECT_NAME)
    
    
    
    binvars=["x","Mh","z","pTtot","xF","Mx"]
    for binvar in binvars:
        plot_channels(data,[#["Fall2018Spring2019_RGA_inbending","precut","piplus_piminus",binvar,"","AZI","standard"],
                            ["Fall2018Spring2019_RGA_inbending","precut","piplus_pi0",binvar,"","AZI","splot_sig"],
                            ["Fall2018_RGA_outbending","precut","piminus_pi0",binvar,"","AZI","splot_sig"]
                            ],
                     #make_title=(True if binvar=="x" else False),
                     drop_edges=(True if binvar!="Mx" else False),
                     out_plot=f"{project_dir}/{project_name}/plots/asym/all_7_with_cuts_{binvar}.pdf")

if __name__ == "__main__":
    
    global io
    
    io = DataIO(PROJECT_NAME="pipi0_paper_RGA_only",
                PROJECT_DIR ="/work/clas12/users/gmat/clas12/clas12_dihadrons/projects"
                output_dir="sideband_dependence")


    # Add the required project_name argument
    parser.add_argument('-p', '--project_name', type=str, required=True, help='the name of the project')

    # Add the optional recreate flag
    parser.add_argument('-c', '--create', action='store_true', help='Create the asymmetry yaml (can take 10+ minutes)')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Access the project_name argument
    project_name = args.project_name

    
    
    if create:
        print("YAY")
    main()
