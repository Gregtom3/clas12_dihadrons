import os
import pickle
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
        os.mkdir(outdir)
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

# Define function parameters
# savedir --> Directory where model is stored
# model_type --> Type of model (gbt/xgb)
# suffix --> One of the headers in utils/subdata.json
def import_model(savedir="",
                 model_type="",
                 suffix=""):
    
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
    dfPars = dfPars.sort_values(by="Importance",ascending=False)
    dfPars.to_csv(f"{savedir}/param_importance_{suffix}.csv",index=False)
    