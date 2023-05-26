from binning import *
import uproot
import numpy as np
from tqdm import tqdm

class DataLoader:
    
    def __init__(self, bin_manager=None):
        self.bin_manager = bin_manager
        self.MCmatch_arr = []
        self.rec_passDihadron_arr = []
        self.gen_passDihadron_arr = []
        
    def read_binned_data_from_mc_file(self, root_file):
        
        if self.bin_manager is None:
            print("Error: Bin manager has not been loaded.")
            return
        
        # Open the ROOT file
        file = uproot.open(root_file)
        # Access the TTree
        tree = file["dihadron_cuts"]
        
        # Indices where REC dihadron passed cuts, and was matched to Monte Carlo dihadron
        self.hit_arr    = (np.array(tree["MCmatch"]) == 1) * (np.array(tree["rec_passDihadron"]) == 1)
        # Indices where the REC dihadron did not pass cuts/was not found (MISS)
        self.miss_arr   = (np.array(tree["rec_passDihadron"])==0)
        # Indices where the REC dihadron passes cuts, but no Monte Carlo dihadron matches (FAKE)
        self.fake_arr   = (np.array(tree["MCmatch"]) == 0) * (np.array(tree["rec_passDihadron"]) == 1)
        
        # Reformat the rect and custom names without their prefix
        rect_names = [rect_name.replace("rec_","").replace("gen_","") for rect_name in self.bin_manager.rect_names]
        custom_names = [custom_name.replace("rec_","").replace("gen_","") for custom_name in self.bin_manager.custom_names]
        
        # Get the branches as numpy arrays
        rect_values = {"rec_"+name: np.array(tree["rec_"+name]) for name in rect_names}
        true_rect_values = {"gen_"+name: np.array(tree["gen_"+name]) for name in rect_names}

        custom_values = {"rec_"+name: np.array(tree["rec_"+name]) for name in custom_names}
        true_custom_values = {"gen_"+name: np.array(tree["gen_"+name]) for name in custom_names}
        
        return rect_values , true_rect_values , custom_values , true_custom_values
    
    def read_binned_data_from_data_file(self, root_files):
        if self.bin_manager is None:
            print("Error: Bin manager has not been loaded.")
            return

        rect_values = {}
        custom_values = {}

        # Reformat the rect and custom names without their prefix
        rect_names = [rect_name.replace("rec_","").replace("gen_","") for rect_name in self.bin_manager.rect_names]
        custom_names = [custom_name.replace("rec_","").replace("gen_","") for custom_name in self.bin_manager.custom_names]
        
        for root_file in tqdm(root_files):
            # Open the ROOT file
            file = uproot.open(root_file)
            # Access the TTree
            tree = file["dihadron_cuts"]

            # Get the branches as numpy arrays
            for name in rect_names:
                if name not in rect_values:
                    rect_values["rec_" + name] = []
                rect_values["rec_" + name].extend(np.array(tree["rec_" + name]))

            for name in custom_names:
                if name not in custom_values:
                    custom_values["rec_" + name] = []
                custom_values["rec_" + name].extend(np.array(tree["rec_" + name]))

        # Convert the lists to numpy arrays
        for name in rect_names:
            rect_values["rec_" + name] = np.array(rect_values["rec_" + name])

        for name in custom_names:
            custom_values["rec_" + name] = np.array(custom_values["rec_" + name])

        return rect_values, custom_values