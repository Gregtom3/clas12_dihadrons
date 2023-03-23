import os
import pickle
import random
import numpy as np
import pandas as pd
from catboost import CatBoostClassifier
from xgboost.sklearn import XGBClassifier
from sklearn.ensemble import RandomForestClassifier

#This function creates directories in the given output directory.
def create_dirs(outdir="",
                nn_type="",
                models=[]):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print("Making directory",outdir)
    
    outdir=outdir+"/"+nn_type
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("Making directory",outdir)
    
    for model in models:
        model_type = model[1]
        subdir=outdir+"/"+model_type
        if not os.path.exists(subdir):
            os.mkdir(subdir)
            os.mkdir(subdir+"/eval_plots")
            print("\tMaking subdirectory",subdir)
            print("\t\tMaking subdirectory",subdir+"/eval_plots")
    print("\n")
    
    return outdir

# This function imports a model from a given directory
# based on the model type and suffix
def import_model(*args):
    if(len(args)==1):
        model_path=args[0]
        # Split model path into parts
        model_parts = model_path.split("/")

        # Get savedir
        savedir = "/".join(model_parts[:-1])

        # Get model_type and suffix
        model_filename = model_parts[-1]
        model_filename_parts = model_filename.split("_")
        model_type = model_filename_parts[0]
        suffix = '_'.join(model_filename_parts[2:])
    elif(len(args)==3):
        savedir=args[0]
        model_type=args[1]
        suffix=args[2]
    
    # Check if model type is GBT
    if(model_type=="gbt"):
        
        # Create a CatBoostClassifier object
        model = CatBoostClassifier()
        
        # Load model from given directory
        model.load_model(f"{savedir}/gbt_model_{suffix}")
    
    # Check if model type is XGBoost
    elif(model_type=="xgb"):
        
        # Load model from given directory
        model = pickle.load(open(f"{savedir}/xgb_model_{suffix}.pkl", "rb"))
        
    # Check if model type is Random Forest
    elif(model_type=="rf"):
        
        # Load model from given directory
        model = pickle.load(open(f"{savedir}/rf_model_{suffix}.pkl", "rb"))
    
    # If model type is unknown
    else:
        
        # Print a message
        print("Unknown model type",model_type,"...returning -1...")
        
        # Return -1
        return -1
    
    # Return model
    return model


#This function makes predictions with a given model
#It takes three parameters:
#model - the model that is used for making the predictions
#model_type - the type of model used for making the predictions (can be either "gbt" or "xgb")
#x - the data to be predicted on

def make_predictions(model=0,
                     model_type="",
                     x=[]):
    
    if(model_type=="gbt" or model_type=="xgb" or model_type=="rf"):
        return np.array(model.predict_proba(x)[:,1],dtype=float)
    else:
        print("Unknown model type",model_type,"...returning -1...")
        return -1
    
# This function saves the feature importances of a given model, such as random forest, gradient boosting tree, or XGBoost, into a csv file. 
# It takes in the model object, the model type (rf, gbt, or xgb), the feature names, the savedir, and a suffix as parameters.
def save_feature_importance(model=0,
                            model_type="",
                            feature_names=[],
                            savedir="",
                            suffix=""):
    # Get feature importances
    if(model_type=="rf"):
        importances = model.feature_importances_
    elif(model_type=="gbt"):
        importances = model.get_feature_importance()
    elif(model_type=="xgb"):
        importances = model.feature_importances_
    
    # Save feature importance to dataframe
    dfPars = pd.DataFrame(data={"Parameter": feature_names,"Importance": importances})
    # Sort by most important
    dfPars = dfPars.sort_values(by="Importance",ascending=False)
    # Save to csv
    dfPars.to_csv(f"{savedir}/param_importance_{suffix}.csv",index=False)
    
    

#Function to get the branch name    
def get_branchname(model_path):
    #Split the file path to get the directory path
    directory_path = os.path.dirname(model_path)
    
    #Get a list of all files in the directory
    files_in_directory = os.listdir(directory_path)
    
    #Loop through all the files in the directory
    for file in files_in_directory:
        #Check if the filename matches the pattern
        if file.startswith("model_params"):
            #Open the file and read the first line
            with open(os.path.join(directory_path, file)) as f:
                first_line = f.readline()
                branchname=first_line.strip()
    
    #Check if the model_path has calo or track
    if("/calo/" in model_path):
        branchname+="_calo"
    elif("/track/" in model_path):
        branchname+="_track"
    
    #Return the branchname
    return branchname

# Generate a random sorting of the data if the user wants to only train on a fraction of the sample
def reindex(pool,percent):
    assert(percent>0 and percent<=1)
    X=np.array(pool[0])
    y=np.array(pool[1])
    N=len(X)
    if(percent==1):
        return [X,y]
    else:
        num_points = int(np.ceil(N*percent)) # Number of data points to keep (n)
        idxs = random.sample(range(len(X)),num_points)  # Get (n) random indecies from 0 to N-1 
        return [ X[idxs] , y[idxs] ]