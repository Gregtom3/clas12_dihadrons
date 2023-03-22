import pandas as pd
import numpy as np
import uproot
import yaml
import json
import os
from sklearn.model_selection import train_test_split

# This function processes a dictionary to determine the keys, original keys, and neighbor parameters of the dictionary. It returns a list of keys, a list of original keys, and a list of neighbor parameters.
# Input: uproot TTree
# Outputs: List of key names, list of modified key names, and index of neighbor
def get_keys(u):
    
    m_g = u["m_g"].array(library="np")[0]
    m_ch = u["m_ch"].array(library="np")[0]
    m_nh = u["m_nh"].array(library="np")[0]
    
    keys = u.keys()
    keys = [key for key in keys if not key in ["m_g","m_ch","m_nh"]]
    
    original_keys=[]
    modified_keys=[]
    nns = []
    
    for key in keys:
        M=-1
        # Since the model can train on multiple neighbor parameters, ensure there
        # is a unique input for each neighbor
        if("_gamma" in key):
            M=m_g
        elif("_ch" in key):
            M=m_ch
        elif("_nh" in key):
            M=m_nh
            
        if(M==-1):
            original_keys.append(key)
            modified_keys.append(key)
            nns.append(-1)
        else:
            for i in range(M):
                original_keys.append(key)
                modified_keys.append(key+f"_{i}")
                nns.append(i)
    
    return original_keys,modified_keys,nns

# Creates dataset for training or evaluation
# This function reads in root files, extracts the TTree and converts it into a dataframe, then splits the data into a training and validation set according to the given split ratio and random seed. The output is four objects, X_train, X_validation, y_train, y_validation containing the training and validation features and labels, respectively.

# Inputs:
#   rootfiles = [file1.root, file2.root, ...] created by photonML.C
#   ttree     = "MLInput" or any other TTree name if needed
#   version   = Determines if data must be split for training
#   split_ratio = Ratio for splitting data for training and testing
#   random_seed = Seed for random number generator
def load_data(rootfiles=[""],
              ttree="MLInput",
              version="train",
              split_ratio = 0.75,
              random_seed = 42):
    assert(version=="train" or version=="evaluate")
    
    if(type(rootfiles)!=list):
        rootfiles=[rootfiles]
        print("WARNING: Need to convert <rootfiles> to list")
    
    df=0
    # Loop over rootfiles
    for ifile,rootfile in enumerate(rootfiles):
        # Open the file in uproot
        u = uproot.open(rootfile)
        u = u[ttree]
        
        # Get dataframe params
        branchnames,keys,nns = get_keys(u)
        
        if(ifile==0): 
            # Create dataframe from first file
            # This dataframe will be appended to for each file
            df = pd.DataFrame(columns=keys)
        
        # Create temporary dataframe for each file
        tmp_df = pd.DataFrame(columns=keys)
        
        # Fill temporary dataframe
        for b,k,n in zip(branchnames,keys,nns):
            if(n==-1):
                tmp_df[k]=u[b].array(library="np")
            else:
                tmp_df[k]=np.array(u[b].array(library="ak")[:,n],dtype=float)
        
        # Concatenate with main dataframe
        df=pd.concat([df,tmp_df], ignore_index=True,axis=0)
    # Create dataset for training/evaluation
    X = df.drop("flag",axis=1)
    if(version=="evaluate"):
        return [X,0],[0,0]
    else:
        y=df["flag"]
        X_train, X_validation, y_train, y_validation = train_test_split(X, y, train_size=split_ratio, random_state=random_seed)
        return [X_train, y_train], [X_validation, y_validation]
    
    
    
#This code opens a yaml file, and creates a list of the parameters stored in the file.
#It does this by looping through each item in the yaml file, and storing the name and type of each item, as well as the parameters associated with it.
def load_params(inyaml=""):
    #Open the yaml file
    with open(inyaml) as file:
        #Load the file into a python dictionary
        yaml_dict = yaml.safe_load(file)
    
    #Create an empty list to store the params
    model_params_list=[]
    #Loop through each item in the dictionary
    for name,model in yaml_dict.items():
        #Create a new empty dictionary to store the parameters
        model_params = dict()
        #Store the type of each item in a variable
        _type = model["type"]
        _type = _type.strip()
        #Loop through the other key-value pairs in the dictionary
        for key, value in model.items():
                #If the key is not 'type', store the key-value pair in the model_params dictionary
                if key != 'type':
                    model_params[key] = value
        #Add the name, type and parameters of each item to the model_params_list
        model_params_list.append([name,_type,model_params])
    
    #Return the list of parameters
    return model_params_list

#This function is used to load the files for the machine learning algorithm.
#It takes 3 parameters: rootdir (the directory of the files), SUBDATA (the subdataset to be used for training) and a boolean test value.
#It then opens the subdata.json file and reads the keys from the JSON file.
#The function then iterates through the files in the rootdir and appends those that end with .root and contain "MC" in the file name to the root_files list. 
#If the SUBDATA parameter is not set to "all", the function checks if the file name contains any of the runs for that SUBDATA and only appends if it does.
#Finally, the function prints the amount of files found and returns the root_files list.
def load_files(rootdir="",
               SUBDATA="MC_inbending",
               test=0):
    
    if(test==1):
        return [rootdir+"/MC_3051_0.root"]
    
    fjs = open ('../../utils/subdata.json', "r")
    JSON = json.loads(fjs.read())
    SUBDATA_KEYS=[key for key in JSON.keys()]
    
    root_files = []

    for file in os.listdir(rootdir):
        if (file.endswith(".root") and "MC" in file):
            foundFile=False
            if(SUBDATA!="all"):
                for RUN in JSON[SUBDATA]["Runs"]:
                    if RUN in file:
                        foundFile=True
            else:
                foundFile=True
            if(foundFile):
                root_files.append(rootdir+"/"+file)
                
    print(len(root_files),"root files found for the ML train/test")
    
    return root_files