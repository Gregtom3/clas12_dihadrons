import uproot
import numpy as np
import ROOT
from tqdm import tqdm
from data_io import *
from binning import *
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors

class SmearingMatrix:

    def __init__(self, root_files, bin_manager):
        
        '''
            Initializes the SmearingMatrix class with a list of root files and binning schemes.
        ''' 
        
        if(type(root_files)!=list):
            root_files=[root_files]
            
        self.root_files = root_files
        self.bin_manager = bin_manager
        self.smearing_matrix = np.zeros((self.bin_manager.total_bins, self.bin_manager.total_bins),dtype=int)
        
        self.dataloader = DataLoader(bin_manager = self.bin_manager)
    
        self.true_bins = np.array([]) # True bin of kinematic
        self.reco_bins = np.array([]) # Reco bin of kinematic
    
    def run(self):
        
        '''
            Main function that reads data from files, fills the smearing matrix
        ''' 
        
        ####
        for root_file in tqdm(self.root_files):
            rect_values , true_rect_values , custom_values , true_custom_values = self.dataloader.read_binned_data_from_mc_file(root_file)
            self.fill_smearing_matrix(rect_values , true_rect_values , custom_values , true_custom_values)
        ####
        
        
    def fill_smearing_matrix(self, rect_values , true_rect_values , custom_values , true_custom_values):

        '''
            Fills the smearing matrix based on the bins identified from true and reco data.
        '''
        
        true_bins = self.bin_manager.get_bin_ids(true_rect_values,true_custom_values)
        reco_bins = self.bin_manager.get_bin_ids(rect_values,custom_values)
        
        for idx,_ in enumerate(true_bins):
            
            # Update the smearing matrix
            true_bin = true_bins[idx]
            reco_bin = reco_bins[idx]
            if(self.dataloader.hit_arr[idx]==1):
                self.smearing_matrix[reco_bin, true_bin] += 1
            elif(self.dataloader.miss_arr[idx]==1):
                reco_bins[idx] = -1
            elif(self.dataloader.fake_arr[idx]==1):
                true_bins[idx] = -1
        
        self.true_bins = np.append(self.true_bins,true_bins)
        self.reco_bins = np.append(self.reco_bins,reco_bins)
        
    def draw_smearing_matrix(self, do_log_scale=False, show_overflow=False):
        plt.figure(dpi=200,figsize=(7,7))
        
        M = np.array(self.smearing_matrix,dtype=float)

        if not show_overflow:
            M = M[1:, 1:]

        # Set zero entries to NaN to make them invisible
        M[M == 0] = np.nan
        

        # Set zero entries to white color
        cmap = matplotlib.colormaps["viridis"]
        cmap.set_bad(color="white")
        
        # Set color scale to log scale if required
        if do_log_scale:
            plt.imshow(M, norm=colors.LogNorm())
        else:
            plt.imshow(M)

        # Set large axis labels and figure title
        plt.xlabel("True Bin",fontsize=12)
        plt.ylabel("Reco Bin",fontsize=12)
        plt.title("Smearing Matrix",fontsize=15)

        # Add color bar with label
        plt.colorbar(label="Counts")

        # Show the plot
        plt.show()

        