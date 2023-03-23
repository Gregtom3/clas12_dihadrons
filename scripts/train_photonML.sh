#!/bin/bash
volatile=/volatile/clas12/users/gmat/clas12analysis.sidis.data
NOW=$( date '+%F_%H_%M_%S' )
farmout=/farm_out/gmat/clas12analysis.sidis.data/clas12_dihadrons/$NOW

# Define the function printred, which prints the argument with the color red
printred () {
  echo -e "\e[31m$1\e[0m"
}

# Define the function printblue, which prints the argument with the color blue
printblue () {
  echo -e "\e[34m$1\e[0m"
}

# Define the function printgreen, which prints the argument with the color green
printgreen () {
  echo -e "\e[32m$1\e[0m"
}

mkdir_green () {
    mkdir -p "$1"
    if [ -n "$2" ]; then
        printf "\t%.s" $(seq 1 "$2")
    fi
    printgreen "mkdir $1"
}


# Assign the variable PWD the value of the current working directory
PWD=`pwd`

# If the current working directory's base name is not "clas12_dihadrons"
if [ $(basename $(pwd)) != "clas12_dihadrons" ]; then
  # Print an error message
  printred "This script must be run from the clas12_dihadrons directory"
  # Exit
  exit
fi

# Prompt the user for a project name
PROJECT_DIR=$PWD/projects
echo "Available projects"
ls $PROJECT_DIR

read -p "Please enter a project name: " PROJECT_NAME

# Check if the project exists in the project directory
if [ -d "$PROJECT_DIR/$PROJECT_NAME" ]; then
  echo "Project exists"
  # Set the project directory to the selected project
  PROJECT_DIR=$PROJECT_DIR/$PROJECT_NAME
  MODEL_DIR=$PROJECT_DIR/models
else
  # Print an error message and exit the program
  echo "Project does not exist...Aborting..."
  exit 1
fi

# Prompt the user for model_params list
PARAMS_DIR=$PWD/machine_learning/photonID/params_folder
echo "Available param lists"
ls $PARAMS_DIR

read -p "Please enter a params list: " PARAMS_NAME

# Check if the file exists
if [ ! -f $PARAMS_DIR/$PARAMS_NAME ]; then
  echo "$PARAMS_DIR/$PARAMS_NAME does not exist. Exiting..."
  exit 1
else
  PARAMS="${PARAMS_DIR}/${PARAMS_NAME}"
fi





VOLATILE_DIR="$volatile/clas12_dihadrons/projects/$PROJECT_NAME"
DATA_DIR="$VOLATILE_DIR/data/raw"
FARMOUT_DIR="$farmout"

# Create folders within FARMOUT_DIR
mkdir_green "$FARMOUT_DIR/log"
mkdir_green "$FARMOUT_DIR/err"
mkdir_green "$FARMOUT_DIR/slurm"

# Create pion pid pairs for the ML portion
pion_pairs=("piplus_pi0" "piminus_pi0" "pi0_pi0")
declare -A pion_pairs_pids
pion_pairs_pids[0,0]=211
pion_pairs_pids[0,1]=111

pion_pairs_pids[1,0]=-211
pion_pairs_pids[1,1]=111

pion_pairs_pids[2,0]=111
pion_pairs_pids[2,1]=111

# Create separate trainable models for the different Monte Carlos
monte_carlos=("MC_RGA_inbending" "MC_RGA_outbending")

# Create separate trainable models using either calo or track training
nn_types=("calo" "track")

# Create the training scripts
for pid_pair in "${pion_pairs[@]}"; do
  for mc in "${monte_carlos[@]}"; do
    for nn_type in "${nn_types[@]}"; do
      rootdir="$DATA_DIR/$pid_pair"
      
      slurmslurm=$FARMOUT_DIR/slurm/train_photonML_${pid_pair}_${mc}_${nn_type}.slurm

      touch -f $slurmslurm

      cat >> $slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition gpu
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=job_photonML_${pid_pair}_${mc}_${nn_type}
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:TitanRTX:1
#SBATCH --time=24:00:00
#SBATCH --output=$FARMOUT_DIR/log/train_photonML_${pid_pair}_${mc}_${nn_type}.out
#SBATCH --error=$FARMOUT_DIR/err/train_photonML_${pid_pair}_${mc}_${nn_type}.err
/u/apps/python3/3.8.7/bin/python3 $PWD/machine_learning/photonID/train.py "$rootdir" "$mc" "$PARAMS" "$nn_type" "$MODEL_DIR/photonID/$pid_pair"
EOF
    echo "Submitting slurm job for ${pid_pair} , ${mc} , ${nn_type}"
    sbatch --quiet $slurmslurm
    done
  done
done


