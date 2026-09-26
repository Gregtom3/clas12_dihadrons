import uproot
import numpy as np
import ROOT
ROOT.gSystem.Load("/home/gmat/link_to_clas12/packages/RooUnfold/libRooUnfold.so")

class AcceptanceMatrix:

    def __init__(self, root_files, rect_bins=None, rect_names=None, custom_bins=None, custom_names=None, method="bayes"):
        
        '''
            Initializes the AcceptanceMatrix class with a list of root files and binning schemes.
        ''' 
        
        if(type(root_files)!=list):
            root_files=[root_files]
            
        self.root_files = root_files
        
        self.rect_bins = self._set_default_rect_bins(rect_bins) # List of rectangular bin schemes
        self.rect_names = self._set_default_rect_names(rect_names) # List of rectangular bin names 
        self.custom_bins = self._set_default_custom_bins(custom_bins) # List of custom bins
        self.custom_names = self._set_default_custom_names(custom_names) # List of variables used in custom bins
        self.total_bins = self._set_total_bins() # Total number of bins, +1 for underflow/overflow
        
        self.acceptance_matrix = np.zeros((self.total_bins, self.total_bins),dtype=int)

        self.true_vector = np.zeros(self.total_bins)
        self.reco_vector = np.zeros(self.total_bins)
    
        self.true_hist = ROOT.TH1D("ht","",self.total_bins-1,1.0,1.0*self.total_bins) # Not counting over/underflow bins
        self.reco_hist = ROOT.TH1D("hr","",self.total_bins-1,1.0,1.0*self.total_bins)
        self.response = ROOT.RooUnfoldResponse (self.total_bins-1,1.0,1.0*self.total_bins)
        self.method = method
        self.hUnfold = None # Unfolded histogram
        
    def _set_total_bins(self):
        # Calculate the total number of bins
        total_bins = int(np.prod([len(edges) - 1 for edges in self.rect_bins]) * len(self.custom_bins)) + 1
        return total_bins
    
    def _set_default_rect_bins(self, rect_bins):
        
        '''
            Sets the default rectangular bins if not provided by the user.
        ''' 
        
        if rect_bins is None:
            self.no_rect_bins = True
            return [[-999, 999]]
        self.no_rect_bins = False
        return rect_bins

    def _set_default_rect_names(self, rect_names):
        
        '''
            Sets the default rectangular names if not provided by the user.
        ''' 
        
        if rect_names is None:
            return ["x"]
        return rect_names

    def _set_default_custom_bins(self, custom_bins):
        
        '''
            Sets the default custom bins if not provided by the user.
        ''' 
        
        if custom_bins is None:
            self.no_custom_bins = True
            return ["([x]>-999)"]
        self.no_custom_bins = False
        return custom_bins

    
    
    def _set_default_custom_names(self, custom_names):
        
        '''
            Sets the default custom names if not provided by the user.
        ''' 
        
        if custom_names is None:
            return ["x"]
        return custom_names
     
        
        
    def read_data_from_file(self, root_file):
        
        '''
            Reads data from the ROOT file and initializes the acceptance matrix.
        ''' 
        
        # Open the ROOT file
        file = uproot.open(root_file)
        # Access the TTree
        tree = file["dihadron_cuts"]
    
        # Get the branches as numpy arrays
        rect_values = {name: np.array(tree[name]) for name in self.rect_names}
        true_rect_values = {name: np.array(tree["true"+name]) for name in self.rect_names}

        custom_values = {name: np.array(tree[name]) for name in self.custom_names}
        true_custom_values = {name: np.array(tree["true"+name]) for name in self.custom_names}
        
        return rect_values , true_rect_values , custom_values , true_custom_values
    
    
    
    
    def run(self):
        
        '''
            Main function that reads data from files, fills the acceptance matrix, and sets vectors.
        ''' 
        
        ####
        for root_file in self.root_files:
            rect_values , true_rect_values , custom_values , true_custom_values = self.read_data_from_file(root_file)
            self.fill_acceptance_matrix(rect_values , true_rect_values , custom_values , true_custom_values)
        ####
        
        self.build_vectors()
        self.build_histos()    
        self.unfold()
     
    
    
    def build_histos(self):
        for i,(x,xt) in enumerate(zip(self.reco_vector, self.true_vector)):
            if i == 0:
                continue
            self.reco_hist.Fill(i,x)
            self.true_hist.Fill(i,xt)
       
    
    
    def unfold(self):
        if self.method == "bayes":
            uf = ROOT.RooUnfoldBayes(self.response, self.reco_hist, 4);    #  OR\
        elif self.method == "svd":
            uf= ROOT.RooUnfoldSvd (self.response, self.reco_hist,20);
            
        self.hUnfold = uf.Hunfold();
        
        
        
    def build_vectors(self):
        
        '''
            Sets the true and reco vectors based on the acceptance matrix.
        '''
        
        self.true_vector = np.sum(self.acceptance_matrix, axis=1)
        self.reco_vector = np.sum(self.acceptance_matrix, axis=0)
        
        
    def fill_acceptance_matrix(self, rect_values , true_rect_values , custom_values , true_custom_values):

        '''
            Fills the acceptance matrix based on the bins identified from true and reco data.
        '''
        
        true_bins = self.get_bin_ids(true_rect_values,true_custom_values)
        reco_bins = self.get_bin_ids(rect_values,custom_values)
        
        for true_bin, reco_bin in zip(true_bins, reco_bins):
            # Update the acceptance matrix
            self.acceptance_matrix[true_bin, reco_bin] += 1
            self.response.Fill(reco_bin,true_bin)
        
        
    def get_bin_ids(self, rect_values, custom_values):
        
        '''
            Returns the bin IDs based on the rectangular and custom values provided.
        '''
    
        # Convert the rect values to just the list of np.arrays
        rect_values = [rect_values[key] for key in rect_values]
        
        # Get the ids for the rectangular bins and custom bins separately
        rect_bin_ids = self.get_rect_bin_ids(rect_values)
        custom_bin_ids = self.get_custom_bin_ids(custom_values)

        # Determine the events where there was overflow/underflow
        bad_indices = ((rect_bin_ids==-1)|(custom_bin_ids==-1))
        
        # Set the unique bin id to be 1 --> N bins (0th bin is reserved for underflow/overflow)
        bin_ids = rect_bin_ids * len(self.custom_bins) + custom_bin_ids + 1
        
        # Set the unique bin id to 0 for underflow/overflow
        bin_ids[bad_indices] = 0
        
        return bin_ids
    
    def get_rect_bin_ids(self,values):
        
        '''
            Returns the rectangular bin IDs based on the values provided.
        '''
        
        # Find the corresponding bin index for each value in parallel
        indices = [np.digitize(val, bins) - 1 for val, bins in zip(values, self.rect_bins)]
        
        # Set the overflow/underflow indices to 0 in parallel
        bad_indices = np.any([(index==-1) | (index==len(bins)-1) for index, bins in zip(indices, self.rect_bins)], axis=0)
        
        # Calculate the bin id based on the indices
        factors = np.cumprod([len(edges) - 1 for edges in self.rect_bins[::-1]])[:-1][::-1]
        # For single binnings, set factors to [1]
        try:
            if np.empty(factors):
                factors = np.array([1])
        except:
            factors = np.append(factors,1)
            
        bin_ids = np.sum([index * factor for index, factor in zip(indices, factors)], axis=0)
        
        bin_ids[bad_indices] = -1


        return bin_ids
        
    def get_custom_bin_ids(self,values):
        
        '''
            Returns the custom bin IDs based on the values provided.
        '''
        
        custom_bins = self.custom_bins
        # Replace variables in the bin expression with corresponding arrays
        for key, value in values.items():
            custom_bins=[cb.replace(f"[{key}]", f"values['{key}']") for cb in custom_bins]

        bin_ids = None
        for idx, cb in enumerate(custom_bins):
            
            mask = eval(cb)
            if idx>0:
                bin_ids += np.where(mask, idx, -1) + 1  # Set bin index as True value
            else:
                bin_ids = np.where(mask, idx, -1)  # Set bin index as True value
        return bin_ids
        
    def get_bins(self, bin_id):
        
        '''
            Returns the bins for a given bin_id.
        '''
        
        if bin_id == 0:
            raise ValueError("bin_id cannot be 0")

        bin_id -= 1  # Offset due to overflow/underflow bins
        rect_bin_id, custom_bin_id = divmod(bin_id, len(self.custom_bins))

        # Get the indices for the rectangular bins
        factors = np.cumprod([len(edges) - 1 for edges in self.rect_bins[::-1]])[:-1][::-1]
        indices = []
        for factor in factors:
            index, rect_bin_id = divmod(rect_bin_id, factor)
            indices.append(index)
        indices.append(rect_bin_id)
        bin_lows = [self.rect_bins[i][indices[i]] for i in range(len(indices))]
        bin_highs = [self.rect_bins[i][indices[i]+1] for i in range(len(indices))]
        rect_bins = [f"([{self.rect_names[i]}]>{bin_lows[i]})&([{self.rect_names[i]}]<{bin_highs[i]})" for i in range(len(indices))]
        #rect_bins = [self.rect_bins[i][indices[i]:indices[i]+2] for i in range(len(indices))]
        
        # Get the bin for the custom bins
        custom_bins = self.custom_bins[custom_bin_id]
    

        if(self.no_rect_bins):
            return custom_bins
        elif(self.no_custom_bins):
            return rect_bins
        return rect_bins, custom_bins


    def convert_to_rect_custom_values(self, values):
        
        '''
            Converts values to rectangular and custom values based on the bin names.
        '''
        
        rect_values = {name: values[name] for name in self.rect_names}
        custom_values = {name: values[name] for name in self.custom_names}
        return rect_values, custom_values
    
    def get_unique_id(self, values):
        
        '''
            Returns a unique bin ID for a given set of values.
        '''
        
        rect_values, custom_values = self.convert_to_rect_custom_values(values)
        return self.get_bin_ids(rect_values,custom_values)[0]
    
    
class CustomBinFactory:
    def __init__(self, pars):
        assert(type(pars) == list)  # Ensure pars is a list
        self.pars = pars  # Store the list of parameters
        self.bins = []  # Initialize an empty list to store custom bins
        self.curves = {}  # Initialize an empty dictionary to store curve definitions
        
    def replace_elements(self, input_string):
        # Replace elements in the input string with brackets
        for element in self.pars:
            input_string = input_string.replace(element, f'[{element}]')
        
        if "[" not in input_string:
            raise ValueError("ERROR: Unable to find variables in custom bin string")
            # If no brackets are found in the modified string, raise an error
        
        return input_string
    
    def add_curve(self, curve_name, s):
        # Add a curve definition to the curves dictionary
        self.curves[curve_name] = self.replace_elements(s)
        # Use 'self.replace_elements' to call the method
    
    def make_bin(self, curve_names, *lines):
        new_bin = []  # Initialize a new bin list
        
        for name in self.curves:
            # Iterate over the defined curves
            if name not in curve_names:
                continue
                # Skip curves not present in curve_names
            
            new_bin.append("(" + self.curves[name] + ")")
            # Add the curve definition to the new bin list
        
        for line in lines:
            # Iterate over the additional lines
            line = self.replace_elements(line)
            # Replace elements in the line
            new_bin.append("(" + line + ")")
            # Add the modified line to the new bin list
        
        new_bin = "&".join(new_bin)
        # Join the new bin list with '&' as the separator
        
        self.bins.append(new_bin)
        # Add the new bin to the bins list
        
    def get_custom_bins(self):
        return self.bins
        # Return the list of custom bins
    
    def get_custom_names(self):
        return self.pars
        # Return the list of parameter names
    
    
            
        