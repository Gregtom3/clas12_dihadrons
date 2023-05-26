from utils import *
import numpy as np


class BinManager:
    def __init__(self):
        
        # Declare a default variable to be a "dummy" bin
        self.default_variable = "rec_x"    
        self.default_range = [-99999,99999] 
        
        self.rect_names = [self.default_variable]
        self.rect_bins = [self.default_range]
        self.custom_names = [self.default_variable]
        self.custom_bins = [f"([{self.default_variable}]>{self.default_range[0]}) & \
                              ([{self.default_variable}]<{self.default_range[1]})"]
        
        self.total_bins = 1
        self.has_rect_bins = False
        self.has_custom_bins = False
        
    def load_factory(self, factory):
        if type(factory)==RectBinFactory:
            self.rect_names = factory.pars
            self.rect_bins  = factory.bins
            self.has_rect_bins = True
            
        elif type(factory)==CustomBinFactory:
            self.custom_names = factory.pars
            self.custom_bins  = factory.bins
            self.has_custom_bins = True
        else:
            raise ValueError("Unknown factory type",type(factory))
        
        self.set_total_bins() # Update number of bins
        
    def set_total_bins(self):
        # Calculate the total number of bins
        self.total_bins = int(np.prod([len(edges) - 1 for edges in self.rect_bins]) * len(self.custom_bins)) + 1
    
    def get_bin_ids(self, rect_values, custom_values):
        
        '''
            Returns the bin IDs based on the rectangular and custom values provided.
            rect_values and custom_values are dicts
            see data_io.py for details
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
        
        for key, _ in values.items():
            if "gen_" in key:
                custom_bins = [cb.replace("rec_", "gen_") for cb in custom_bins]
            break
            
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
    
    def get_all_rect_bin_ids_from_rect_bin_id(self,rect_bin_id):
        # Get the indices for the rectangular bins
        factors = np.cumprod([len(edges) - 1 for edges in self.rect_bins[::-1]])[:-1][::-1]
        indices = []
        for factor in factors:
            index, rect_bin_id = divmod(rect_bin_id, factor)
            indices.append(index)
        indices.append(rect_bin_id)
        return indices
    
    def get_bin_ids_from_unique_id(self,bin_id):
        
        if bin_id == 0:
            raise ValueError("bin_id cannot be 0")

        bin_id -= 1  # Offset due to overflow/underflow bins
        rect_bin_id, custom_bin_id = divmod(bin_id, len(self.custom_bins))
        return rect_bin_id,custom_bin_id
    
    def get_bins_from_unique_id(self, bin_id):
        
        rect_bin_id, custom_bin_id = self.get_bin_ids_from_unique_id(bin_id)

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

        # Get the bin for the custom bins
        custom_bins = self.custom_bins[custom_bin_id]
    
        return rect_bins, custom_bins

    def get_tcut_from_unique_id(self,unique_id):
        
        '''
            Given a bin's unique id, return the TCut for drawing from within that bin
        '''
        
        rect_bins, custom_bin = self.get_bins_from_unique_id(unique_id)
        rect_bins = "&".join(rect_bins)
        rect_bin = rect_bins.replace("[","").replace("]","").replace("&","&&")
        custom_bin = custom_bin.replace("[","").replace("]","").replace("&","&&")
        tcut = "&&".join([rect_bin,custom_bin])
        return tcut

    def convert_to_rect_custom_values(self, values):
        
        '''
            Converts values to rectangular and custom values based on the bin names.
        '''
        
        rect_values = {name: values[name] for name in self.rect_names}
        custom_values = {name: values[name] for name in self.custom_names}
        return rect_values, custom_values
    
    def get_unique_id_from_values(self, values):
        
        '''
            Returns a unique bin ID for a given set of values.
            Need to be in the form of a dict
        '''
        
        rect_values, custom_values = self.convert_to_rect_custom_values(values)
        return self.get_bin_ids(rect_values,custom_values)[0]
    
    def print_bins_from_values(self,values):
        unique_id = self.get_unique_id_from_values(values)
        s = ''.join([key+f"={value} , " for key, value in values.items()])
        print("Unique ID for",s,"--->",unique_id)
        self.print_bins_from_unique_id(unique_id)
        
    def print_bins_from_unique_id(self,unique_id):
        bin_edges = self.get_bins_from_unique_id(unique_id)
        print("Bin edges for Unique id",unique_id,"\n"+25*"==")
        self.print_bins_from_bin_edges(bin_edges)
       
    def print_bins_from_bin_edges(self,bin_edges):
        for be in bin_edges:
            if(type(be)==list):
                for _be in be:
                    print(_be)
            else:
                print(be,"\n")

class RectBinFactory:
    def __init__(self):
        self.pars = []
        self.bins = []  # Initialize an empty list to store custom bins
    
    def make_bins(self, par, bins):
        if par in self.pars:
            raise ValueError("Cannot use the same parameter twice")
        
        self.pars.append(par)
        self.bins.append( np.array(bins, dtype=float) ) 



class CustomBinFactory:
    def __init__(self, pars):
        assert(type(pars) == list)  # Ensure pars is a list
        self.pars = pars  # Store the list of parameters
        self.bins = []  # Initialize an empty list to store custom bins
        self.curves = {}  # Initialize an empty dictionary to store curve definitions
        
    
    def add_curve(self, curve_name, s):
        # Add a curve definition to the curves dictionary
        self.curves[curve_name] = replace_elements_with_brackets(s, self.pars)
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
            line = replace_elements_with_brackets(line,self.pars)
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

        
        
    