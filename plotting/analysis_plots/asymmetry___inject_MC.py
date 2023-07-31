import argparse
from tools___inject import *
from tools___io import *
from tools___etc import *
# Create parser
parser = argparse.ArgumentParser(
                prog='asymmetry__inject_MC.py',
                description='Inject asymmetries into dihadron Monte Carlo and extract them (primarily using sideband technique)')

parser.add_argument('--infile', type=str, help='Monte Carlo file for injection')
parser.add_argument('--project_name' , type=str, help='Name of the Analysis Project', required = True)
parser.add_argument('--n_trials', type=int, help='Number of injection trials (default: 10)', default=10)
parser.add_argument('--n_cpus', type=int, help='Max number of CPUs running at a time (default: 4)', default=4)
parser.add_argument('--program', type=argparse.FileType('r'), help='Input file path designating injection setup', default = "./injection_files/default.txt")

args = parser.parse_args()

#####################################################################
# Program: asymmetry___inject_MC.py
# Author:  Gregory Matousek
# Date:    7/25/2023
#
# Purpose: Inject asymmetries into dihadron Monte Carlo and extract them
# Requirements: Completed data processing pipeline to obtain merged TTrees
#####################################################################

            
            
def run_project_trial(i):
    np.random.seed(project.seeds[i])
    project.create_injector()
    yaml_file = f"{project.project_loc}/trial_{i}.yaml"
    project.injector.run(print_to_yaml=True,yaml_name=yaml_file)
    print("Saving to YAML -->","/".join(yaml_file.split["/"][-5:]))
    
def main():
    
    infile = args.infile
    pion_pair = pull_pair_from_string(infile) # ex: "piplus_pi0"
    
    global project
    # Create injection project with multiple trials + parallel processing
    project = InjectionProject(#project_title  = args.inject_project_name, # Injection project name
                               project_title = args.program.name.split("/")[-1]+"_"+pion_pair,
                               infile        = infile,
                               project_loc   = io.out, # Output directory
                               pion_pair     = pion_pair,
                               n_trials      = args.n_trials, # Total number of trials
                               n_cpus        = args.n_cpus,
                               injector_program_file = args.program)   # Max number of CPUs
    
    with ProcessPoolExecutor(max_workers=project.n_cpus) as executor:
        executor.map(run_project_trial,  range(project.n_trials))

    project.make_plots()
    
    
if __name__ == "__main__":
    
    global io
    
    io = DataIO(PROJECT_NAME="pipi0_paper_RGA_only",
                file_names=[args.infile],
                output_dir=f"inject_use_true_phis/{args.program.name.split('/')[-1]+'_'+pull_pair_from_string(args.infile)}")
    
    main()
