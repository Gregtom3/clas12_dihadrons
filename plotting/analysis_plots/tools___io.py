import os
import ROOT
import matplotlib.pyplot as plt

class DataIO:
    def __init__(self, PROJECT_NAME="", PROJECT_DIR="", file_paths=[], file_names=[], ttree_names=[], output_dir=""):
        self.PROJECT_NAME = PROJECT_NAME
        self.PROJECT_DIR = PROJECT_DIR
        self.file_paths = file_paths
        self.file_names = file_names
        self.ttree_names = ttree_names
        self.out = os.getcwd()+"/output_plots/"+self.PROJECT_NAME+"/"+output_dir
        if not os.path.exists(self.out):
            os.makedirs(self.out)
            print("Created new directory:",self.out)

    # When called, the program will create multiple pion pair subdirectories within output_dir
    def create_multiple_pion_pair_dirs(self):
        for pion_pair in ["piplus_piplus", "piplus_pi0", "piplus_piminus", "pi0_pi0", "piminus_pi0", "piminus_piminus"]:
            out = self.out+"/"+pion_pair
            if not os.path.exists(out):
                os.makedirs(out)
                print("Created new directory:",out)
    

    
        
    