#Import packages 
import pandas as pd 
import numpy as np 
from sklearn.ensemble import RandomForestClassifier
import pickle
import os

# Load data into train and test set
def train(train_pool,
          validation_pool,
          params, # Dictionary of model parameters
          outdir="",
          suffix=""):
    
    X_train = np.array(train_pool[0],dtype=float)
    y_train = np.array(train_pool[1],dtype=float)
    
    # Create the model
    model = RandomForestClassifier(**params)
    
    # Fit the model
    model.fit(X_train, y_train)
    
    # Save the model using pickle
    pickle.dump(model,open(outdir+f"/rf_model_{suffix}.pkl","wb"))
    
if __name__ == "__main__":
    train()