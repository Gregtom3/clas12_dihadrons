#!/bin/bash

PROJECT_NAME="$1"
nFiles=$2
nEvents=$3
model_list="$4"
ml_branch="$5" # <model name>_<track/calo>
                       # See model_list in machine_learning/photonID/params_folder
configuration="$6"





hl="-------------------------------------------------------------------"
function wait_for_jobs() {
    local job_name=$1
    local jobsLeft=-999
    
    while [ $jobsLeft -ne 0 ]
    do
        read jobsLeft <<< $(echo "$(squeue -u gmat --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R")" | grep $job_name | awk 'END{print NR}')
        echo "Jobs remaining: " $jobsLeft
        sleep 30
    done
}

# Define a function that calls another function with user inputs
function call_function {
  local function_name=$1 # Set the function name
  local wait_str=$2 # Wait string for jobs
  local inputs=("${@:3}") # Get the user inputs as an array
  
  # Call the function with the inputs
  echo $hl
  echo "${function_name}" "${inputs[@]}"
  echo $hl
  sleep 2
  bash "${function_name}" "${inputs[@]}"

  # Wait until the slurm is finished
  wait_for_jobs $wait_str
}

#
module unload root
module load clas12/pro
#

# 1. create_project.sh
# (Converts hipo files to TTrees with cuts)
# ----------------------------------------
in1=("-o" $PROJECT_NAME $nFiles $nEvents $configuration)
wait1="job_hipo2tree"
func1="./scripts/create_project.sh"

# Call the function
call_function "${func1}" "${wait1}" "${in1[@]}"

# 2. preproess_photonML.sh
# (Creates TTree for machine learning)
# ----------------------------------------
in2=($PROJECT_NAME)
wait2="job_photonML_"
func2="./scripts/preprocess_photonML.sh"

# Call the function
call_function "${func2}" "${wait2}" "${in2[@]}"

# 3. train_photonML.sh
# (Performs training on Monte Carlo)
# ----------------------------------------
in3=($PROJECT_NAME $model_list $configuration)
wait3="job_photonMLtrain_"
func3="./scripts/train_photonML.sh"

# Call the function
call_function "${func3}" "${wait3}" "${in3[@]}"



# 4. predict_photonML.sh
# (Uses trained model to make predictions on MC and Data)
# ----------------------------------------
in4=($PROJECT_NAME $ml_branch $configuration)
wait4="job_photonMLpredict_"
func4="./scripts/predict_photonML.sh"

# Call the function
call_function "${func4}" "${wait4}" "${in4[@]}"




# 5. form_dihadrons.sh
# (Reads TTrees with ML to form dihadron pairs)
# ----------------------------------------
in5=($PROJECT_NAME $ml_branch)
wait5="job_dihadron_"
func5="./scripts/form_dihadrons.sh"

# Call the function
call_function "${func5}" "${wait5}" "${in5[@]}"


# 6. merge_dihadrons.sh
# (Merge the TTrees by data versions)
# ----------------------------------------
in6=($PROJECT_NAME $configuration)
wait6="job_merge_"
func6="./scripts/merge_dihadrons.sh"

# Call the function
call_function "${func6}" "${wait6}" "${in6[@]}"

