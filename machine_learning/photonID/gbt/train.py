import pandas as pd
import numpy as np 
from sklearn.model_selection import train_test_split
from catboost import CatBoostClassifier, Pool, metrics, cv
from sklearn.metrics import accuracy_score
import catboost
import matplotlib.pyplot as plt

# Load data into train and test set
def train(train_pool,
          validation_pool,
          params, # Dictionary of model parameters
          outdir=".",
          suffix=""):
    
    
    X_train = train_pool[0]
    y_train = train_pool[1]
    X_validation = validation_pool[0]
    y_validation = validation_pool[1]
    
    numeric_train_pool = Pool(X_train, y_train)
    numeric_val_pool = Pool(X_validation, y_validation)

    model = CatBoostClassifier(**params,
                            custom_loss=[metrics.Accuracy()], 
                            random_seed=42,
                            task_type="GPU",
                            devices='0:1')
    model.fit(numeric_train_pool, verbose=1,eval_set=numeric_val_pool)

    model.save_model(outdir+f"/gbt_model_{suffix}")

if __name__ == "__main__":
    train()


