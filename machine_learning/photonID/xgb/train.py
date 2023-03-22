import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn.metrics import roc_curve, auc, confusion_matrix
import matplotlib.pyplot as plt
import pickle
import os

# Load data into train and test set
def train(train_pool,
          validation_pool,
          params, # Dictionary of model parameters
          outdir="",
          suffix=""):
    
    X_train = train_pool[0]
    y_train = train_pool[1]
    X_validation = validation_pool[0]
    y_validation = validation_pool[1]
    
    # Create the model
    model = XGBClassifier(**params)
    
    # Fit the model
    model.fit(X_train, y_train)
    
    # Save the model using pickle
    pickle.dump(model,open(outdir+f"/xgb_model_{suffix}.pkl","wb"))
    
if __name__ == "__main__":
    train()