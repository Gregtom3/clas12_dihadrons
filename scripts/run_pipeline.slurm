#!/bin/bash

PROJECT_NAME="ana_test"
nFiles=5
nEvents=10000
model_list="model_params.yaml"
ml_branch="gbt_0_calo" # <model name>_<track/calo>
                       # See model_list in machine_learning/photonID/params_folder


#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=job_pipeline_${PROJECT_NAME}
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=out_${PROJECT_NAME}.out
#SBATCH --error=err_${PROJECT_NAME}.err

bash ./scripts/pipeline.sh ${PROJECT_NAME} ${nFiles} ${nEvents} ${model_list} ${ml_branch}