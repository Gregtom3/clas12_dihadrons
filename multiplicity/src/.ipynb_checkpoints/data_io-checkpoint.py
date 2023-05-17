from binning import *
import uproot
import numpy as np

class DataLoader:
    
    def __init__(self, bin_manager=None):
        self.bin_manager = bin_manager
        
    def read_binned_data_from_file(self, root_file):
        
        if self.bin_manager is None:
            print("Error: Bin manager has not been loaded.")
            return
        
        # Open the ROOT file
        file = uproot.open(root_file)
        # Access the TTree
        tree = file["dihadron_cuts"]
        
        # Get the indices where the full event matched with Monte Carlo
        good_idx    = np.array(tree["MCmatch"]) == 1
        
        # Get the branches as numpy arrays
        rect_values = {name: np.array(tree[name])[good_idx] for name in self.bin_manager.rect_names}
        true_rect_values = {name: np.array(tree["true"+name])[good_idx] for name in self.bin_manager.rect_names}

        custom_values = {name: np.array(tree[name])[good_idx] for name in self.bin_manager.custom_names}
        true_custom_values = {name: np.array(tree["true"+name])[good_idx] for name in self.bin_manager.custom_names}
        
       
        return rect_values , true_rect_values , custom_values , true_custom_values