import sys
import ROOT
import numpy as np
import array
import os
from tqdm import tqdm
from tools.dataloader import *
from tools.handler import *
from catboost import CatBoostClassifier
from xgboost.sklearn import XGBClassifier
from sklearn.ensemble import RandomForestClassifier

def predict(rootdir="/volatile/clas12/users/gmat/clas12analysis.sidis.data/clas12_dihadrons/projects/ana_v0/data/pi0_pi0",
            SUBDATA="Fall2018_RGA_inbending",
            model_path="/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/ana_test/models/photonID/piplus_pi0/calo/gbt_0"):
    
    # Load in rootfiles for analysis
    rootfiles = load_files(rootdir=rootdir,
                           SUBDATA=SUBDATA)
    
    # Import the trained model
    model = import_model(model_path)
    
    # Determine the branch name for the classifier output
    branchname = get_branchname(model_path)
    print("Branchname -->",branchname)
    
    # Determine the ttree containing the MLInput based on the model
    if("/calo/" in model_path):
        ttree="MLInput_calo"
    elif("/track/" in model_path):
        ttree="MLInput_track"
    
    # For loop over each file
    for ifile,rootfile in enumerate(tqdm(rootfiles)):
        # Load MLInput data
        X=load_data(rootfiles=[rootfile],
                    ttree=ttree,
                    version="predict")

        # No data was successfully loaded
        # We return -1 when data loading fails
        if(type(X)==int):
            continue

        # Make signal predictions
        prob=model.predict_proba(X)[:,1]
        
        # Load EventTree
        # Here we will create a new weights branch so that the photon-per-photon classification can be stored
        tfile = ROOT.TFile(rootfile, "UPDATE")
        tree = tfile.Get("EventTree")

        # create a new branch for the tree
        weights = array.array('d',100*[0.0])
        weight_branch=tree.Branch(branchname, weights,branchname+'[Nmax]/D')
        
        # Loop over the events in the EventTree and set the weights array accordingly
        # If the particle is not a photon, set weight to 1
        # If the particle is a photon, set the weight from the prediction
        k=0 # An incrementing variable to step to the next element in "prob" after a photon is found
        N=tree.GetEntries()
        for i in range(N):
            tree.GetEntry(i)
            pid=np.array(tree.pid)
            for j,PID in enumerate(pid):
                if(PID==22):
                    weights[j]=prob[k]
                    k+=1
                else:
                    weights[j]=1
            weight_branch.Fill() # Fill the branch
        
        # Write the TTree and close the TFile
        tree.Write("EventTree",1) # the "1" forces an overwrite of the previous ttree
        tfile.Close()
        
if __name__ == "__main__":
    args = [arg for arg in sys.argv[1:]]
    predict(*args)    
        
        
            
