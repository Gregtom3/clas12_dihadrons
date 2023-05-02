import sys
import os
from tools.dataloader import *
from tools.handler import *
from tools.plotter import *
from random_forest.train import train as rf_train
from gbt.train import train as gbt_train
from xgb.train import train as xgb_train

# ML imports
# from catboost import CatBoostClassifier, Pool, metrics, cv
# from sklearn.metrics import accuracy_score
# from xgboost.sklearn import XGBClassifier
# from sklearn.metrics import roc_curve, auc, confusion_matrix


def train(rootdir = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/clas12_dihadrons/projects/ana_vrgc/data/piplus_pi0",
          SUBDATA= "MC_RGC",
          yamlfile = "/work/clas12/users/gmat/clas12/clas12_dihadrons/machine_learning/photonID/params_folder/model_params_gbt_only_full.yaml",
          nn_type  = "calo", # calo or track (use either calorimeter info or track info to determine nearest neighbors)
          outdir   = "/work/clas12/users/gmat/clas12/clas12_dihadrons/projects/ana_vrgc/models/photonID/piplus_pi0",
          input_model = None):
    
    # Load the parameters from the yamlfile
    if(input_model==None):
        models = load_params(yamlfile)
    else:
        models = [input_model]

    # Create the output directories to store the models
    # Sets outdir=outdir/nn_type
    outdir=create_dirs(outdir = outdir,
                       nn_type = nn_type,
                       models = models)
    
    # Load the rootfiles for the learning
    # SUBDATA --> see utils/subdata.json
    #         --> Must be one of the headers (ex: MC_RGA_inbending)
    rootfiles = load_files(rootdir=rootdir,
                           SUBDATA=SUBDATA)
    print(len(rootfiles),"found for training")
    
    # Split the data into a training and validation set
    if(nn_type=="calo"):
        ttree="MLInput_calo"
    elif(nn_type=="track"):
        ttree="MLInput_track"
    else:
        print("Unknown nn_type",nn_type,"...defaulting to MLInput_calo...")
        ttree="MLInput_calo"
    
    # Load in data
    print("Loading data for",len(rootfiles),"root files")
    train_pool, validation_pool = load_data(rootfiles = rootfiles,
                                            ttree     = ttree,
                                            version="train",
                                            split_ratio = 0.75,
                                            random_seed = 42)
    
    # For each of the models in the yamlfile, perform the training
    print("Starting training")
    for model in models:
        model_name = model[0]
        model_type = model[1]
        model_percentage = model[2]
        model_params = model[3]
        
        # Define savedir
        savedir=outdir+"/"+model_name
        
        print("\n\nProcessing model",model_name,"with type",model_type,"\n\n")
        
        # Based on the model percentage, only train on a fraction of the data
        tpool = reindex(train_pool,model_percentage)
        vpool = reindex(validation_pool,model_percentage)
        
        # Determine which architecture to train with
        # The model is saved to "savedir" which can be imported later
        if(model_type=="gbt"):
            gbt_train(tpool,
                      vpool,
                      model_params,
                      savedir,
                      SUBDATA)
        
        elif(model_type=="xgb"):
            xgb_train(tpool,
                      vpool,
                      model_params,
                      savedir,
                      SUBDATA)
        
        elif(model_type=="rf"):
            rf_train(tpool,
                     vpool,
                     model_params,
                     savedir,
                     SUBDATA)
        else:
            print("Unknown model type:",model_type,"...skipping...")
            continue
        
        # Write the model parameters to a text file in the dictionary
        with open(savedir+"/model_params.txt", "w") as f:
            f.write(model_name+"\n")
            for key, value in model_params.items():
                f.write(key + ": " + str(value) + "\n")

        # Import the model and perform predictions with the validation set
        trained_model = import_model(savedir,model_type, SUBDATA)
        
        # Make predictions from the model
        X_validation = validation_pool[0]
        y_validation = validation_pool[1]
        predictions = make_predictions(trained_model,model_type,X_validation)
        
        # Create plots and save them to the model directories
        make_plots(X_validation,y_validation,predictions,savedir,SUBDATA)
        
        # Save the parameter importances
        feature_names=X_validation.keys()
        save_feature_importance(trained_model,model_type,feature_names,savedir,SUBDATA)

if __name__ == "__main__":
    args = [arg for arg in sys.argv[1:]]
    train(*args)